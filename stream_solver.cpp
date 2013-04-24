#include "stream_solver.h"

stream_solver::stream_solver(const MatrixXd &r, const MatrixXd &z, const MatrixXd &t,
                             const int &bn, const int &sn, const int &stn,
                             const double &in_p0, const double &in_t0, const double &out_p, const MatrixXd &cir,
                             const VectorXd &circulation_in, const VectorXd &entropy_in,
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
            ||tmp_i_2!=cir.cols()
            ||tmp_i_1!=circulation_in.size()
            ||cir(0,0)!=circulation_in(0)
            ||cir(tmp_i_1-1,0)!=circulation_in(tmp_i_1-1)
            ||tmp_i_1!=entropy_in.size())
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
    inlet_total_pressure=in_p0;
    inlet_total_temperature=in_t0;
    outlet_pressure=out_p;
    Cp=heat_capacity_ratio*gas_constant/(heat_capacity_ratio-1);
    total_enthalpy_inlet.resize(stream_number);
    total_enthalpy_inlet.fill(Cp*inlet_total_temperature);
    circulation_inlet=circulation_in;
    entropy_inlet=entropy_in;
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
    meridian_stream_length.col(0).fill(0);
    for(int n1=0;n1<stream_number;++n1)
    {
        //in radial turbine, z increase with flow, and r reduce with flow
        //assume flow direction of turbine wheel
        //at leading edge, flow direction is oposite of r Vector2d(0,-1)
        //at trailing edge, flow direction is along z Vector2d(1,0)
        spline s(z_axial.row(n1).transpose(),radius.row(n1).transpose(),
                 spline_endpoint_direction_2d(Vector2d(0,-1)),
                 spline_endpoint_direction_2d(Vector2d(1,0)));
        tmp_d_1=s.get_dif_1();
        tmp_d_2=s.get_dif_2();
        tmp_d_3=s.get_length();
        meridian_stream_direction_z.row(n1)=tmp_d_1.col(0).transpose();
        meridian_stream_direction_r.row(n1)=tmp_d_1.col(1).transpose();
        tmp_arr_d_1=tmp_d_1.col(0).array();
        tmp_arr_d_2=tmp_d_1.col(1).array();
        tmp_arr_d_3=tmp_d_2.col(0).array();
        tmp_arr_d_4=tmp_d_2.col(1).array();
        meridian_stream_curvature.row(n1)=(tmp_arr_d_4/tmp_arr_d_3
                                           *std::sqrt(std::pow((tmp_arr_d_1*tmp_arr_d_1/(tmp_arr_d_1*tmp_arr_d_1+tmp_arr_d_2*tmp_arr_d_2)),3))).transpose();
        meridian_stream_length.row(n1)=tmp_d_3.transpose();
    }
}

void stream_solver::calculate_area()
{
    MatrixXd arc_length, tmp1, tmp2;
    tmp1=radius.topRows(stream_number-1)-radius.bottomRows(stream_number-1);
    tmp2=z_axial.topRows(stream_number-1)-z_axial.bottomRows(stream_number-1);
    delta_q.row(0).fill(0);
    delta_q.bottomRows(stream_number-1)=(tmp1.cwiseProduct(tmp1)+tmp2.cwiseProduct(tmp2)).cwiseSqrt();
    arc_length=radius.cwiseProduct(MatrixXd::Constant(stream_number,station_number,2*M_PI)-blade_number*thickness);
    meridian_area=(arc_length.topRows(stream_number-1)+arc_length.bottomRows(stream_number-1)).cwiseProduct(delta_q.bottomRows(stream_number-1))/2;
    q=delta_q.colwise().sum().colwise().replicate(stream_number);
}

void stream_solver::flow_field_initialization()
{
    radius=radius_original;
    z_axial=z_axial_original;
    circulation=circulation_original;
    thickness=thickness_original;
    calculate_stream_directions();
    calculate_area();
    beta.resizeLike(radius);
    beta.fill(0);
    theta.resizeLike(radius);
    theta.fill(M_PI/6);
    pressure.resizeLike(radius);
    total_pressure.resizeLike(radius);
    temperature.resizeLike(radius);
    density.resizeLike(radius);
    enthalpy.resizeLike(radius);
    total_enthalpy.resizeLike(radius);
    entropy.resizeLike(radius);
    relative_speed_r.resizeLike(radius);
    relative_speed_z.resizeLike(radius);
    relative_speed_theta.resizeLike(radius);
    double in_p, in_t, in_den;
    thermaldynamic_equations::total_to_static(inlet_total_pressure, inlet_total_temperature,
                                              meridian_area.col(0).sum()*cos(M_PI/6),
                                              mass_flow_rate,gas_constant,heat_capacity_ratio,
                                              in_p,in_t,in_den);
    pressure.fill(in_p);
    total_pressure.fill(inlet_total_pressure);
    temperature.fill(in_t);
    density.fill(in_den);
    enthalpy.fill(Cp*in_t);
    total_enthalpy.fill(Cp*inlet_total_temperature);
    entropy.col(0)=entropy_inlet;
    meridian_stream_direction_r.resizeLike(radius);
    meridian_stream_direction_z.resizeLike(radius);
    meridian_stream_curvature.resizeLike(radius);
    meridian_stream_length.resize(stream_number,station_number-1);
    efficiency_grid.resizeLike(meridian_stream_length);
    meridian_area.resize(stream_number-1,station_number);
    dtheta_dm.resizeLike(radius);
    dcirculation_dm.resizeLike(radius);
    dwm_dm.resizeLike(radius);
    //initialize speed field
    MatrixXd tmp_mat_1, tmp_mat_2, spd_mat;
    spd_mat=(MatrixXd::Constant(1,station_number,mass_flow_rate).cwiseQuotient(meridian_area.colwise().sum())/density(0,0)).colwise().replicate(stream_number);
    tmp_mat_1=(meridian_stream_direction_r.cwiseProduct(meridian_stream_direction_r)
               +meridian_stream_direction_z.cwiseProduct(meridian_stream_direction_z)).cwiseSqrt();
    tmp_mat_1=(z_axial.row(stream_number-1)-z_axial.row(0)).colwise().replicate(stream_number);
    tmp_mat_2=(radius.row(stream_number-1)-radius.row(0)).colwise().replicate(stream_number);
    relative_speed_m=spd_mat.cwiseQuotient(math_ext::cos_subtract(meridian_stream_direction_z, meridian_stream_direction_r,
                                                                  tmp_mat_2,tmp_mat_1));
    relative_speed_r=relative_speed_m.cwiseProduct(meridian_stream_direction_r).cwiseQuotient(tmp_mat_1);
    relative_speed_z=relative_speed_m.cwiseProduct(meridian_stream_direction_z).cwiseQuotient(tmp_mat_1);
    rothalpy.col(0)=MatrixXd::Constant(stream_number,1,Cp*inlet_total_temperature)-circulation.col(0)*rotate_speed;
}

void stream_solver::calculate_s2m()
{
    curvature_centrifugal_force=calculate_curvature_centrifugal();
    interpolate_circulation();
    calculate_dtheta_dm();
    calculate_theta();
    pressure_grandiant_force=calculate_pressure_gradiant();
    thermal_grandiant_force=calculate_pressure_gradiant();
    //page 252, equation 3-77
    MatrixXd dwm_dq=curvature_centrifugal_force.cwiseProduct(relative_speed_m)
            +pressure_grandiant_force+thermal_grandiant_force.cwiseQuotient(relative_speed_m);
    calculate_efficiency_grid();
    calculate_thermaldynamic();
    for(int n1=1;n1<station_number;++n1)
    {
        double tmp_mf, residence, wm_hub=relative_speed_m(0,0);
        do
        {
            tmp_mf=calculate_station_mass_flow(wm_hub,dwm_dq.col(n1),n1);
            residence=(tmp_mf-mass_flow_rate)/mass_flow_rate;
            wm_hub*=1+residence;
        }while(std::fabs(residence)<0.001);
    }
}

MatrixXd stream_solver::calculate_curvature_centrifugal()//coefficient A, page 252, equation 3-77a
{
    MatrixXd cos_a_substract_b, tmp1, tmp2;
    tmp1=(z_axial.row(stream_number-1)-z_axial.row(0)).colwise().replicate(stream_number);
    tmp2=(radius.row(stream_number-1)-radius.row(0)).colwise().replicate(stream_number);
    cos_a_substract_b=math_ext::cos_subtract(meridian_stream_direction_z, meridian_stream_direction_r,
                                             tmp2,tmp1);
    return(cos_a_substract_b.cwiseProduct(meridian_stream_curvature));
}

void stream_solver::interpolate_circulation()
{
    //because I use quasi-orthogonal mesh for computation
    //and I move node along quasi-orthogonal line(movement is constrained)
    //so interpolate along quasi-orthogonal line only
    MatrixXd tmp(stream_number,station_number);
    delta_circulation=(circulation_original.row(stream_number-1)-circulation_original.row(0)).colwise().replicate(stream_number);
    tmp.row(0).fill(0);
    for(int n1=1;n1<stream_number;++n1)
    {
        tmp.row(n1)=tmp.row(n1-1)+delta_q.row(n1);
    }
    circulation=circulation_original.row(0).colwise().replicate(stream_number)
            +tmp.cwiseProduct(delta_circulation).cwiseQuotient(q);
}

void stream_solver::calculate_theta()//page 252, equation 3-77d
{
    MatrixXd tmp1=(dtheta_dm.leftCols(station_number-1)+dtheta_dm.rightCols(station_number-1))/2;
    delta_theta=tmp1.cwiseProduct(meridian_stream_length);
    for(int n1=1;n1<station_number;++n1)
    {
        theta.col(n1)=theta.col(n1-1)+delta_theta.col(n1-1);
    }
}

void stream_solver::calculate_dtheta_dm()
{
    MatrixXd tmp1, tmp2;
    tmp1=rotate_speed*radius.cwiseProduct(radius);
    tmp2=relative_speed_m.cwiseProduct(radius).cwiseProduct(radius);
    dtheta_dm=(circulation-tmp1).cwiseQuotient(tmp2);
}

MatrixXd stream_solver::calculate_pressure_gradiant()//coefficient B, page 252, equation 3-77b
{
    MatrixXd tmp1, tmp2(stream_number,station_number);
    VectorXd tmp_vec_1(station_number), tmp_vec_2(station_number);
    tmp_vec_1(0)=0;
    //stream line length m is sorted by spanwise
    //so calculate dcirculation_dm and dwm_dm by spanwise
    for(int n1=0;n1<stream_number;++n1)
    {
        for(int n2=1;n2<station_number;++n2)
        {
            tmp_vec_1(n2)=tmp_vec_1(n2-1)+meridian_stream_lenth(n1,n2-1);
        }
        spline s(tmp_vec_1,circulation.row(n1));
        tmp1=s.get_dif_1();
        dcirculation_dm.row(n1)=tmp1.col(1).cwiseQuotient(tmp1.col(0)).transpose();
        tmp_vec_2=(relative_speed_r.cwiseProduct(relative_speed_r)
                   +relative_speed_z.cwiseProduct(relative_speed_z)
                   ).cwiseSqrt();
        s.input(tmp_vec_1,tmp_vec_2);
        tmp1=s.get_dif_1();
        dwm_dm.row(n1)=tmp1.col(1).cwiseQuotient(tmp1.col(0)).transpose();
    }
    //spanwise length q is sorted by flow direction
    //so calculate dtheta_dq and dcirculation_dq by stations
    tmp2.row(0).fill(0);
    for(int n1=1;n1<stream_number;++n1)
    {
        tmp2.row(n1)=tmp2.row(n1-1)+delta_q.row(n1);
    }
    for(int n1=0;n1<station_number;++n1)
    {
        spline s(tmp2.col(n1),theta.col(n1));
        tmp1=s.get_dif_1();
        dtheta_dq.col(n1)=tmp1.col(1).cwiseQuotient(tmp1.col(0));
        s.input(tmp2.col(n1),circulation.col(n1));
        tmp1=s.get_dif_1();
        dcirculation_dq.col(n1)=tmp1.col(1).cwiseQuotient(tmp1.col(0));
    }
    tmp1=(z_axial.row(stream_number-1)-z_axial.row(0)).colwise().replicate(stream_number);
    tmp2=(radius.row(stream_number-1)-radius.row(0)).colwise().replicate(stream_number);
    return(dcirculation_dm.cwiseProduct(dtheta_dq)
           -dcirculation_dq.cwiseProduct(dtheta_dm)
           +dwm_dm.cwiseProduct(math_ext::sin_subtract(meridian_stream_direction_z, meridian_stream_direction_r,
                                                       tmp2,tmp1)));
}

MatrixXd stream_solver::calculate_thermal_gradiant()//coefficient C, page 252, equation 3-77c
{
    MatrixXd tmp1(stream_number,station_number), tmp2,
            total_enthalpy_gradient(stream_number,station_number),
            circulation_gradient(stream_number,station_number),
            entropy_gradient(stream_number,station_number);
    for(int n1=1;n1<stream_number;++n1)
    {
        tmp1.row(n1)=tmp1.row(n1-1)+delta_q.row(n1);
    }
    for(int n1=0;n1<station_number;++n1)
    {
        spline s(tmp1.col(n1),total_enthalpy_inlet);
        tmp2=s.get_dif_1();
        total_enthalpy_gradient.col(n1)=tmp2.rowwise().norm();
        s.input(tmp1.col(n1),circulation_inlet);
        tmp2=s.get_dif_1();
        circulation_gradient.col(n1)=tmp2.rowwise().norm();
        s.input(tmp1.col(n1),entropy_inlet);
        tmp2=s.get_dif_1();
        entropy_gradient.col(n1)=tmp2.rowwise().norm();
    }
    circulation_gradient*=rotate_speed;
    entropy_gradient=entropy_gradient.cwiseProduct(temperature);
    return(total_enthalpy_gradient+circulation_gradient+entropy_gradient);
}

double stream_solver::calculate_station_mass_flow(const double &wm0, const VectorXd &dif_1, const int &station)
//page 229, equation 3-41
{
    //page 231, equation 3-43 to 3-45
    VectorXd wm=math_ext::runge_kutta(wm0,dif_1,delta_q.block(1,station,stream_number-1,1)),
            spd_square(wm.size()), tmp1, tmp2, tmp_p, tmp_t, tmp_rho;
    const int n=dif_1.size();
    //by the assumed wm0, we can get all wm at the station
    //by the wm at the station and circulation, we can get speed at the station
    //we got total pressure and total temperature in calculate_thermaldynamic()
    //we can get static pressure, temperature and density by total parameter and speed
    //finally we can get mass flow by the wm and density
    //in page 264, the 3rd block is this function, but it's not explained detailed in the literature
    tmp1=circulation.col(station).cwiseQuotient(radius.col(station));
    spd_square=wm.cwiseProduct(wm)+tmp1.cwiseProduct(tmp1);
    tmp_t=total_temperature.col(station)-0.5*spd_square/Cp;
    tmp_p=total_pressure.col(station).cwiseQuotient(
                MatrixXd::Ones(n,1)+spd_square.cwiseQuotient(tmp_t)/(2*gas_constant));
    tmp_rho=tmp_p.cwiseQuotient(tmp_t)/gas_constant;
    tmp1=wm.cwiseProduct(tmp_rho);
    tmp2=(tmp1.head(n-1)+tmp1.tail(n-1))/2;
    return(tmp2.cwiseProduct(meridian_area.col(station)).sum());
}

void stream_solver::calculate_efficiency_grid()
{
    //Jiang Zikang, page 359, (6)
    VectorXd stream_length=meridian_stream_length.rowwise().sum(),
            tmp_vec_1=MatrixXd::Ones(stream_number,1),
            tmp_vec_2=MatrixXd::Constant(stream_number,1,wheel_efficiency);
    for(int n1=0;n1<station_number-1;++n1)
    {
        efficiency_grid.col(n1)=MatrixXd::Ones(stream_number,1)
                -(tmp_vec_1-tmp_vec_2).cwiseQuotient(stream_length).cwiseProduct(meridian_stream_length.col(n1));
    }
}

void stream_solver::calculate_thermaldynamic()
{
    //Jiang Zikang, page 358, equation 14-94
    for(int n1=1;n1<station_number;++n1)
    {
        total_enthalpy.col(n1)=total_enthalpy.col(n1-1)+rotate_speed
                *(circulation.col(n1)-circulation.col(n1-1));
    }
    //Jiang Zikang, page 358, equation 14-96
    total_temperature=total_enthalpy/Cp;
    ArrayXd tmp_arr_1(stream_number,1);
    for(int n1=1;n1<station_number;++n1)
    {
        //Jiang Zikang, page 359, equation 14-98
        tmp_arr_1=(ArrayXd::Ones(stream_number,1)+efficiency_grid.col(n1-1).cwiseProduct(
                       total_temperature.col(n1).cwiseQuotient(total_temperature.col(n1-1))
                       -MatrixXd::Ones(stream_number,1)
                       ).array()
                   ).pow(heat_capacity_ratio/(heat_capacity_ratio-1));
        total_pressure.col(n1)=total_pressure.col(n1-1).cwiseProduct(tmp_arr_1.matrix());
        //Jiang Zikang, page 359, equation 14-99
        tmp_arr_1=total_temperature.col(n1).cwiseQuotient(total_temperature.col(n1-1)).array().log()*Cp
                -total_pressure.col(n1).cwiseQuotient(total_pressure.col(n1-1)).array().log()*gas_constant;
        entropy.col(n1)=entropy.col(n1-1)+tmp_arr_1.matrix();
    }

}
