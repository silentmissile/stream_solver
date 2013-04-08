#include "stream_solver.h"

stream_solver::stream_solver(const MatrixXd &r, const MatrixXd &z, const MatrixXd &t,
                             const int &bn, const int &sn, const int &stn,
                             const double &in_p, const double &in_t, const double &out_p, const MatrixXd &cir,
                             const double &mf, const double &rs, const MatrixXd &eff,
                             const double &R, const double &gamma)
{
    //check input data
    int tmp_i_1, tmp_i_2;
    tmp_i_1=r.rows();
    tmp_i_2=r.cols();
    if(tmp_i_1!=z.rows()
            ||tmp_i_2!=z.cols()
            ||tmp_i_1!=t.rows()
            ||tmp_i_2!=t.cols()
            ||tmp_i_1!=eff.rows()
            ||tmp_i_2!=eff.cols()
            ||tmp_i_1!=cir.rows()
            ||tmp_i_2!=cir.cols())
        return;
    //import input data
    stream_number=tmp_i_1;
    station_number=tmp_i_2;
    radius_original=r;
    z_axial_original=z;
    thickness_original=t;
    circulation_original=cir;
    blade_number=bn;
    stream_number=sn;
    station_number=stn;
    inlet_pressure=in_p;
    inlet_temperature=in_t;
    outlet_pressure=out_p;
    mass_flow_rate=mf;
    rotate_speed=rs;
    wheel_efficiency=eff;
    gas_constant=R;
    heat_capacity_ratio=gamma;
    //initialize flow field
    flow_field_initialization();
}

void stream_solver::calculate_stream_directions()
{
    MatrixXd tmp_d_1(station_number,2), tmp_d_2(station_number,2);
    VectorXd tmp_d_3;
    ArrayXd tmp_arr_d_1, tmp_arr_d_2, tmp_arr_d_3, tmp_arr_d_4;
    for(int n1=0;n1<stream_number;++n1)
    {
        //in radial turbine, z increase with flow, and r reduce with flow
        //assume flow direction of turbine wheel
        //at leading edge, flow direction is oposite of r Vector2d(0,-1)
        //at trailing edge, flow direction is along z Vector2d(1,0)
        if(!spline::spline_parametric(z_axial.row(n1),radius.row(n1),
                                      Vector2d(0,-1),Vector2d(1,0),
                                      tmp_d_1,tmp_d_2,tmp_d_3))
            return;
        meridian_stream_direction_z.row(n1)=tmp_d_1.col(0).transpose();
        meridian_stream_direction_r.row(n1)=tmp_d_1.col(1).transpose();
        tmp_arr_d_1=tmp_d_1.col(0).array();
        tmp_arr_d_2=tmp_d_1.col(1).array();
        tmp_arr_d_3=tmp_d_2.col(0).array();
        tmp_arr_d_4=tmp_d_2.col(1).array();
        meridian_stream_curvature.row(n1)=(tmp_arr_d_4/tmp_arr_d_3
                                           *std::sqrt(std::pow((tmp_arr_d_1*tmp_arr_d_1/(tmp_arr_d_1*tmp_arr_d_1+tmp_arr_d_2*tmp_arr_d_2)),3))).transpose();
        meridian_stream_lenth.row(n1)=tmp_d_3.transpose();
    }
}

void stream_solver::calculate_area()
{
    ArrayXd od_distance, arc_length;
    ArrayXd tmp_arr_1, tmp_arr_2;
    tmp_arr_1=radius.topRows(stream_number-1)-radius.bottomRows(stream_number-1);
    tmp_arr_2=z_axial.topRows(stream_number-1)-z_axial.bottomRows(stream_number-1);
    od_distance=std::sqrt(tmp_arr_1*tmp_arr_1+tmp_arr_2*tmp_arr_2);
    arc_length=radius.cwiseProduct(MatrixXd::Constant(stream_number,station_number,2*M_PI)-blade_number*thickness);
    meridian_area=(arc_length.topRows(stream_number-1)+arc_length.bottomRows(stream_number-1))/2*od_distance;
}

void stream_solver::flow_field_initialization()
{
    radius=radius_original;
    z_axial=z_axial_original;
    circulation=circulation_original;
    thickness=thickness_original;
    beta.resizeLike(radius);
    beta.fill(0);
    theta.resizeLike(radius);
    theta.fill(0);
    pressure.resizeLike(radius);
    temperature.resizeLike(radius);
    density.resizeLike(radius);
    enthalpy.resizeLike(radius);
    entropy.resizeLike(radius);
    relative_speed_r.resizeLike(radius);
    relative_speed_z.resizeLike(radius);
    relative_speed_theta.resizeLike(radius);
    pressure.col(0).setConstant(inlet_pressure);
    temperature.col(0).setConstant(inlet_temperature);
    density.col(0).setConstant(inlet_pressure/gas_constant/inlet_temperature);
    entropy.col(0).setConstant(0);
    meridian_stream_direction_r.resizeLike(radius);
    meridian_stream_direction_z.resizeLike(radius);
    meridian_stream_curvature.resizeLike(radius);
    meridian_stream_lenth.resize(stream_number,station_number-1);
    meridian_area.resize(stream_number-1,station_number);
    calculate_stream_directions();
    calculate_area();
    double rho;
    //initialize speed field
    MatrixXd tmp_mat_1, spd_mat;
    rho=inlet_pressure/inlet_temperature/gas_constant;
    density.fill(rho);
    spd_mat=(MatrixXd::Constant(1,station_number,mass_flow_rate).cwiseQuotient(meridian_area.colwise().sum())/rho).colwise().replicate(stream_number);
    tmp_mat_1=(meridian_stream_direction_r.cwiseProduct(meridian_stream_direction_r)
               +meridian_stream_direction_z.cwiseProduct(meridian_stream_direction_z)).cwiseSqrt();
    relative_speed_r=spd_mat.cwiseProduct(meridian_stream_direction_r).cwiseQuotient(tmp_mat_1);
    relative_speed_z=spd_mat.cwiseProduct(meridian_stream_direction_z).cwiseQuotient(tmp_mat_1);
}

void stream_solver::calculate_s2m()
{
    MatrixXd curvature_centrifugal_force;
    calculate_curvature_centrifugal(curvature_centrifugal_force);
    interpolate_circulation();
}

void stream_solver::calculate_curvature_centrifugal(MatrixXd &res)
{
    MatrixXd cos_a, cos_b, sin_a, sin_b, tmp, tmp1, tmp2;
    tmp=(meridian_stream_direction_r.cwiseProduct(meridian_stream_direction_r)
         +meridian_stream_direction_z.cwiseProduct(meridian_stream_direction_z)).cwiseSqrt();
    cos_a=meridian_stream_direction_r.cwiseQuotient(tmp);
    sin_a=meridian_stream_direction_z.cwiseQuotient(tmp);
    tmp1=z_axial.row(stream_number-1)-z_axial.row(0);
    tmp2=radius.row(stream_number-1)-radius.row(0);
    tmp=(tmp1.cwiseProduct(tmp1)+tmp2.cwiseProduct(tmp2)).cwiseSqrt();
    cos_b=tmp1.cwiseQuotient(tmp).rowwise().replicate(stream_number);
    sin_b=tmp2.cwiseQuotient(tmp).rowwise().replicate(stream_number);
    res=(cos_a.cwiseProduct(cos_b)+sin_a.cwiseProduct(sin_b)).cwiseProduct(meridian_stream_curvature);
}

void stream_solver::interpolate_circulation()
{
    MatrixXd q, delta_q, delta_circulation;
    MatrixXd tmp1, tmp2;
    tmp1=z_axial_original.row(stream_number-1)-z_axial_original.row(0);
    tmp2=radius_original.row(stream_number-1)-radius.row(0);
    q=(tmp1.cwiseProduct(tmp1)+tmp2.cwiseProduct(tmp2)).cwiseSqrt().rowwise().replicate(stream_number);
    tmp1=z_axial-z_axial.row(0).rowwise().replicate(stream_number);
    tmp2=radius-radius.row(0).rowwise().replicate(stream_number);
    delta_q=(tmp1.cwiseProduct(tmp1)+tmp2.cwiseProduct(tmp2)).cwiseSqrt();
    delta_circulation=(circulation_original.row(stream_number-1)-circulation_original.row(0)).rowwise().replicate(stream_number);
    circulation=circulation_original.row(0).rowwise().replicate(stream_number)
            +delta_q.cwiseProduct(delta_circulation).cwiseQuotient(q);
}
