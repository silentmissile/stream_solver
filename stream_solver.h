#ifndef STREAM_SOLVER_H
#define STREAM_SOLVER_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "thermaldynamic_equations.h"
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
                           const double &in_p, const double &in_t, const double &out_p, const MatrixXd &cir,
                           const VectorXd &circulation_in, const VectorXd &entropy_in,
                           const double &mf, const double &rs, const double &eff,
                           const double &R, const double &gamma);
private:
    void flow_field_initialization();
    void calculate_stream_directions();
    void calculate_area();
    void calculate_s2m();
    MatrixXd calculate_curvature_centrifugal();//coefficient A
    MatrixXd calculate_pressure_gradiant();//coefficient B
    MatrixXd calculate_thermal_gradiant();//coefficiant C
    void interpolate_circulation();
    void calculate_theta();
    void calculate_dtheta_dm();
    void calculate_efficiency_grid();
    void calculate_thermaldynamic();
    void interpolate_mass_flow(const int &station, const VectorXd &mf_pre);
    double calculate_station_mass_flow(const double &wm_hub, const VectorXd &dif_1, const int &station, VectorXd &mf_distribution);
    int blade_number, stream_number, station_number;
    //in following geometry parameter matrixes
    //each row is for one stream line, from leading edge to trailing edge
    //in each column, the points is from hub to tip
    //thickness is angular thickness not dimensional thickness
    MatrixXd radius_original, z_axial_original, circulation_original, thickness_original,
    radius, z_axial, thickness, beta, theta,
    meridian_stream_direction_z, meridian_stream_direction_r,
    meridian_stream_curvature, meridian_stream_length,
    meridian_area;
    //boundary conditions
    double inlet_total_pressure, inlet_total_temperature, outlet_pressure, mass_flow_rate;
    MatrixXd circulation;
    double rotate_speed;
    double wheel_efficiency;
    VectorXd total_enthalpy_inlet, circulation_inlet, entropy_inlet, mass_flow_rate_inlet;
    //material property
    double gas_constant, heat_capacity_ratio, Cp;
    //flow parameter for grid
    MatrixXd pressure, total_pressure, temperature, total_temperature, density,
    enthalpy, total_enthalpy, entropy, rothalpy, efficiency_grid,
    relative_speed_r, relative_speed_z, relative_speed_m, relative_speed_theta,
    delta_theta, q, delta_q, delta_circulation,
    curvature_centrifugal_force, pressure_grandiant_force, thermal_grandiant_force,
    dtheta_dm, dcirculation_dm, dwm_dm, dtheta_dq, dcirculation_dq;
};

#endif // STREAM_SOLVER_H
