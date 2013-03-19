#ifndef SPLINE_H
#define SPLINE_H
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;

class spline
{
public:
    spline();
    //spline cubic interpolation
    //for x, y 2D spline
    static bool spline_parametric(const VectorXd &x, const VectorXd &y, VectorXd &dif_1);
    static bool spline_parametric(const VectorXd &x, const VectorXd &y, MatrixXd &dif_1);
    static bool spline_parametric(const VectorXd &x, const VectorXd &y,
                                  MatrixXd &dif_1, MatrixXd &dif_2);
    static bool spline_parametric(const VectorXd &x, const VectorXd &y,
                                  const Vector2d &start_direction, const Vector2d &end_direction,
                                  MatrixXd &dif_1);
    static bool spline_parametric(const VectorXd &x, const VectorXd &y,
                                  const Vector2d &start_direction, const Vector2d &end_direction,
                                  MatrixXd &dif_1, MatrixXd &dif_2);
    static bool spline_parametric(const VectorXd &x, const VectorXd &y,
                                  const Vector2d &start_direction, const Vector2d &end_direction,
                                  MatrixXd &dif_1, MatrixXd &dif_2, VectorXd &len);
    //for x, y, z 3D spline
    static bool spline_parametric(const VectorXd &x, const VectorXd &y, const VectorXd &z,
                                  MatrixXd &dif_1);
    static bool spline_parametric(const VectorXd &x, const VectorXd &y, const VectorXd &z,
                                  MatrixXd &dif_1, MatrixXd &dif_2);
    static bool spline_parametric(const VectorXd &x, const VectorXd &y, const VectorXd &z,
                                  const Vector3d &start_direction, const Vector3d &end_direction,
                                  MatrixXd &dif_1);
    static bool spline_parametric(const VectorXd &x, const VectorXd &y, const VectorXd &z,
                                  const Vector3d &start_direction, const Vector3d &end_direction,
                                  MatrixXd &dif_1, MatrixXd &dif_2);
    static bool spline_parametric(const VectorXd &x, const VectorXd &y, const VectorXd &z,
                                  const Vector3d &start_direction, const Vector3d &end_direction,
                                  MatrixXd &dif_1, MatrixXd &dif_2, VectorXd &len);
private:
    //spline core arithmatic
    static bool spline_arithmatic(const VectorXd &h, const VectorXd &y,
                                  VectorXd &dif_1);
    static bool spline_arithmatic(const VectorXd &h, const VectorXd &y,
                                  const double &start_dif_1, const double &end_dif_1,
                                  VectorXd &dif_1);
    static bool spline_arithmatic(const VectorXd &h, const VectorXd &y,
                                  VectorXd &dif_1, VectorXd &dif_2);
    static bool spline_arithmatic(const VectorXd &h, const VectorXd &y,
                                  const double &start_dif_1, const double &end_dif_1,
                                  VectorXd &dif_1, VectorXd &dif_2);
    static bool thomas_solver3(const VectorXd &a, const VectorXd &b, const VectorXd &d, VectorXd &c);
};

#endif // SPLINE_H
