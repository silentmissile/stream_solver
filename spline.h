#ifndef SPLINE_H
#define SPLINE_H
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;

class spline_endpoint
{
public:
    enum order{first, second};
    spline_endpoint(order d, double v);
    order dif;
    double value;
};

class spline_endpoint_direction_2d
{
public:
    spline_endpoint_direction_2d(const Vector2d &in);
    spline_endpoint_direction_2d(const spline_endpoint &x_dir, const spline_endpoint &y_dir);
    spline_endpoint x_direction, y_direction;
};

class spline_endpoint_direction_3d
{
public:
    spline_endpoint_direction_3d(const Vector3d &in);
    spline_endpoint_direction_3d(const spline_endpoint &x_dir,
                                 const spline_endpoint &y_dir,
                                 const spline_endpoint &z_dir);
    spline_endpoint x_direction, y_direction, z_direction;
};

class spline
{
public:
    enum dim{two, three};
    explicit spline(const VectorXd &x_in, const VectorXd &y_in,
                    const spline_endpoint_direction_2d &start_direction
                    =spline_endpoint_direction_2d(spline_endpoint(spline_endpoint::second,0.0),
                                                  spline_endpoint(spline_endpoint::second,0.0)),
                    const spline_endpoint_direction_2d &end_direction
                    =spline_endpoint_direction_2d(spline_endpoint(spline_endpoint::second,0.0),
                                                  spline_endpoint(spline_endpoint::second,0.0)));
    explicit spline(const VectorXd &x_in, const VectorXd &y_in, const VectorXd &z_in,
                    const spline_endpoint_direction_3d &start_direction
                    =spline_endpoint_direction_3d(spline_endpoint(spline_endpoint::second,0.0),
                                                  spline_endpoint(spline_endpoint::second,0.0),
                                                  spline_endpoint(spline_endpoint::second,0.0)),
                    const spline_endpoint_direction_3d &end_direction
                    =spline_endpoint_direction_3d(spline_endpoint(spline_endpoint::second,0.0),
                                                  spline_endpoint(spline_endpoint::second,0.0),
                                                  spline_endpoint(spline_endpoint::second,0.0)));
    void input(const VectorXd &x_in, const VectorXd &y_in,
               const spline_endpoint_direction_2d &start_direction
               =spline_endpoint_direction_2d(spline_endpoint(spline_endpoint::second,0.0),
                                             spline_endpoint(spline_endpoint::second,0.0)),
               const spline_endpoint_direction_2d &end_direction
               =spline_endpoint_direction_2d(spline_endpoint(spline_endpoint::second,0.0),
                                             spline_endpoint(spline_endpoint::second,0.0)));
    void input(const VectorXd &x_in, const VectorXd &y_in, const VectorXd &z_in,
               const spline_endpoint_direction_3d &start_direction
               =spline_endpoint_direction_3d(spline_endpoint(spline_endpoint::second,0.0),
                                             spline_endpoint(spline_endpoint::second,0.0),
                                             spline_endpoint(spline_endpoint::second,0.0)),
               const spline_endpoint_direction_3d &end_direction
               =spline_endpoint_direction_3d(spline_endpoint(spline_endpoint::second,0.0),
                                             spline_endpoint(spline_endpoint::second,0.0),
                                             spline_endpoint(spline_endpoint::second,0.0)));
    inline MatrixXd get_dif_2() const{return(dif_2);}
    MatrixXd get_dif_1();
    MatrixXd get_coefficient();
    VectorXd get_length();
private:
    dim dimension;
    VectorXd x, y, z;
    VectorXd h, length;
    MatrixXd dif_1, dif_2, coefficient;
    void calculate_coefficient();
    //spline core arithmatic
    static bool spline_coefficient_arithmatic(const VectorXd &h, const VectorXd &y, const VectorXd &dif_2,
                                               MatrixXd &coeff);
    static bool spline_arithmatic(const VectorXd &h, const VectorXd &y,
                                  const spline_endpoint &start_point, const spline_endpoint &end_point,
                                  VectorXd &dif_2);
    static bool thomas_solver3(const VectorXd &a, const VectorXd &b, const VectorXd &d, VectorXd &c);
    static VectorXd calculate_chord(const VectorXd &x, const VectorXd &y);
    static VectorXd calculate_chord(const VectorXd &x, const VectorXd &y, const VectorXd &z);
};

#endif // SPLINE_H
