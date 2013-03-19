#include "stream_solver.h"

stream_solver::stream_solver(const MatrixXd &r, const MatrixXd &z, const MatrixXd &t,
                             const int &bn, const int &sn, const int &stn,
                             const double &in_p, const double &in_t, const double &in_spd, const double &out_p,
                             const double &rs, const MatrixXd &eff, const double &R, const double &gamma)
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
            ||tmp_i_2!=eff.cols())
        return;
    //import input data
    stream_number=tmp_i_1;
    station_number=tmp_i_2;
    radius=r;
    z_axial=z;
    thickness=t;
    blade_number=bn;
    stream_number=sn;
    station_number=stn;
    inlet_pressure=in_p;
    inlet_temperature=in_t;
    outlet_pressure=out_p;
    mass_flow_rate=in_spd;
    rotate_speed=rs;
    wheel_efficiency=eff;
    gas_constant=R;
    heat_capacity_ratio=gamma;
    //initialize flow field
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
    VectorXd tmp_vec_1, tmp_vec_2;
    tmp_vec_1=(thickness.col(0).head(stream_number-1)+thickness.col(0).head(stream_number-1))/2;
    tmp_vec_2=z_axial.col(0).tail(stream_number-1)+z_axial.col(0).head(stream_number-1);
    double area=tmp_vec_2.cwiseProduct(MatrixXd::Constant(stream_number-1,1,2*M_PI)-tmp_vec_1*blade_number).sum()*radius(0,0);
    //in radial turbine theory, optimized inlet angle is between 20deg to 40deg
    //here I set it at middle: M_PI/6
    mass_flow_rate=density(0,0)*inlet_speed*cos(M_PI/6)*area;
    //calculate stream directions
    meridian_stream_direction_r.resizeLike(radius);
    meridian_stream_direction_z.resizeLike(radius);
    meridian_stream_curvature.resizeLike(radius);
    meridian_stream_lenth.resize(stream_number,station_number-1);
    stream_direction();
}

void stream_solver::stream_direction()
{
    MatrixXd tmp_d_1(station_number,2), tmp_d_2(station_number,2);
    VectorXd tmp_d_3;
    ArrayXd tmp_arr_d_1, tmp_arr_d_2, tmp_arr_d_3, tmp_arr_d_4;
    for(int n1=0;n1<stream_number;++n1)
    {
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
