#include "math_ext.h"

math_ext::math_ext()
{
}

MatrixXd math_ext::cos_plus(const MatrixXd &dx1, const MatrixXd &dy1,
                            const MatrixXd &dx2, const MatrixXd &dy2)
{
    MatrixXd cos_1, cos_2, sin_1, sin_2;
    sin_and_cos(dx1, dy1, dx2, dy2, sin_1, cos_1, sin_2, cos_2);
    return(cos_1.cwiseProduct(cos_2)-sin_1.cwiseProduct(sin_2));
}

MatrixXd math_ext::cos_subtract(const MatrixXd &dx1, const MatrixXd &dy1,
                                const MatrixXd &dx2, const MatrixXd &dy2)
{
    MatrixXd cos_1, cos_2, sin_1, sin_2;
    sin_and_cos(dx1, dy1, dx2, dy2, sin_1, cos_1, sin_2, cos_2);
    return(cos_1.cwiseProduct(cos_2)+sin_1.cwiseProduct(sin_2));
}

MatrixXd math_ext::sin_plus(const MatrixXd &dx1, const MatrixXd &dy1,
                            const MatrixXd &dx2, const MatrixXd &dy2)
{
    MatrixXd cos_1, cos_2, sin_1, sin_2;
    sin_and_cos(dx1, dy1, dx2, dy2, sin_1, cos_1, sin_2, cos_2);
    return(sin_1.cwiseProduct(cos_2)+cos_1.cwiseProduct(sin_2));
}

MatrixXd math_ext::sin_subtract(const MatrixXd &dx1, const MatrixXd &dy1,
                                const MatrixXd &dx2, const MatrixXd &dy2)
{
    MatrixXd cos_1, cos_2, sin_1, sin_2;
    sin_and_cos(dx1, dy1, dx2, dy2, sin_1, cos_1, sin_2, cos_2);
    return(sin_1.cwiseProduct(cos_2)-cos_1.cwiseProduct(sin_2));
}

void math_ext::sin_and_cos(const MatrixXd &dx1, const MatrixXd &dy1,
                           const MatrixXd &dx2, const MatrixXd &dy2,
                           MatrixXd &sin_1, MatrixXd &cos_1,
                           MatrixXd &sin_2, MatrixXd &cos_2)
{
    MatrixXd tmp;
    tmp=(dy1.cwiseProduct(dy1)
         +dx1.cwiseProduct(dx1)).cwiseSqrt();
    cos_1=dx1.cwiseQuotient(tmp);
    sin_1=dy1.cwiseQuotient(tmp);
    tmp=(dy2.cwiseProduct(dy2)
         +dx2.cwiseProduct(dx2)).cwiseSqrt();
    cos_2=dx2.cwiseQuotient(tmp);
    sin_2=dy2.cwiseQuotient(tmp);
}
