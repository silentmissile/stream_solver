#ifndef STREAM_SOLVER_H
#define STREAM_SOLVER_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include "spline.h"
#include "math_ext.h"
#include "stream_solver_global.h"

using namespace std;
using namespace Eigen;

class STREAM_SOLVERSHARED_EXPORT stream_solver
{
public:
    //stream_solver();
    explicit stream_solver(const MatrixXd &r, const MatrixXd &z, const MatrixXd &t,
                  const int &bn, const int &sn, const int &stn,
                  const double &in_p, const double &in_t, const double &out_p,  const MatrixXd &cir,
                  const double &mf, const double &rs, const MatrixXd &eff,
                  const double &R, const double &gamma);
private:
    void flow_field_initialization();
    void calculate_stream_directions();
    void calculate_area();
    void calculate_s2m();
    MatrixXd calculate_curvature_centrifugal();//coefficient A
    MatrixXd calculate_pressure_grandiant();//coefficient B
    MatrixXd calculate_thermal_grandiant();//coefficiant C
    void interpolate_circulation();
    void calculate_theta();
    void calculate_dtheta_dm();
    int blade_number, stream_number, station_number;
    //in following geometry parameter matrixes
    //each row is for one stream line, from leading edge to trailing edge
    //in each column, the points is from hub to tip
    //thickness is angular thickness not dimensional thickness
    MatrixXd radius_original, z_axial_original, circulation_original, thickness_original;
    MatrixXd radius, z_axial, thickness, beta, theta;
    MatrixXd meridian_stream_direction_z, meridian_stream_direction_r;
    MatrixXd meridian_stream_curvature, meridian_stream_lenth;
    MatrixXd meridian_area;
    //boundary conditions
    double inlet_pressure, inlet_temperature, outlet_pressure, mass_flow_rate;
    MatrixXd circulation;
    double rotate_speed;
    MatrixXd wheel_efficiency;
    //material property
    double gas_constant, heat_capacity_ratio;
    //flow parameter for grid
    MatrixXd pressure, temperature, density, enthalpy, entropy;
    MatrixXd relative_speed_r, relative_speed_z, relative_speed_m, relative_speed_theta;
    MatrixXd delta_theta, q, delta_q, delta_circulation;
    MatrixXd curvature_centrifugal_force, pressure_grandiant_force, thermal_grandiant_force;
    MatrixXd dtheta_dm, dcirculation_dm, dwm_dm, dtheta_dq, dcirculation_dq;
};

#endif // STREAM_SOLVER_H
