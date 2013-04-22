#ifndef MATH_EXT_H
#define MATH_EXT_H
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;

class math_ext
{
public:
    math_ext();
    static MatrixXd cos_plus(const MatrixXd &dx1, const MatrixXd &dy1,
                                   const MatrixXd &dx2, const MatrixXd &dy2);
    static MatrixXd cos_subtract(const MatrixXd &dx1, const MatrixXd &dy1,
                                       const MatrixXd &dx2, const MatrixXd &dy2);
    static MatrixXd sin_plus(const MatrixXd &dx1, const MatrixXd &dy1,
                                   const MatrixXd &dx2, const MatrixXd &dy2);
    static MatrixXd sin_subtract(const MatrixXd &dx1, const MatrixXd &dy1,
                                       const MatrixXd &dx2, const MatrixXd &dy2);
    static VectorXd runge_kutta(const double &x0, const VectorXd &dif_1, const VectorXd &q);
private:
    static void sin_and_cos(const MatrixXd &dx1, const MatrixXd &dy1,
                            const MatrixXd &dx2, const MatrixXd &dy2,
                            MatrixXd &sin_1, MatrixXd &cos_1,
                            MatrixXd &sin_2, MatrixXd &cos_2);
};

#endif // MATH_EXT_H
