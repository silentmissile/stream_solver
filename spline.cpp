#ifndef INTEGRATION_CONSTANT
#define INTEGRATION_CONSTANT 20
#include "spline.h"

spline::spline(const VectorXd &x_in, const VectorXd &y_in,
               const spline_endpoint_direction_2d &start_direction,
               const spline_endpoint_direction_2d &end_direction)
{
    input(x_in, y_in, start_direction,end_direction);
}

spline::spline(const VectorXd &x_in, const VectorXd &y_in, const VectorXd &z_in,
               const spline_endpoint_direction_3d &start_direction,
               const spline_endpoint_direction_3d &end_direction)
{
    input(x_in, y_in, z_in, start_direction,end_direction);
}

void spline::input(const VectorXd &x_in, const VectorXd &y_in,
               const spline_endpoint_direction_2d &start_direction,
               const spline_endpoint_direction_2d &end_direction)
{
    const int n=x_in.size();
    if(n!=y_in.size())
        return;
    dimension=spline::two;
    x=x_in;
    y=y_in;
    h=calculate_chord(x,y);
    VectorXd tmp0, tmp1;
    if(!spline_arithmatic(h,x,start_direction.x_direction,end_direction.x_direction,tmp0)
            ||!spline_arithmatic(h,y,start_direction.y_direction,end_direction.y_direction,tmp1))
        return;
    dif_2.resize(n,2);
    dif_2.col(0)=tmp0;
    dif_2.col(1)=tmp1;
}

void spline::input(const VectorXd &x_in, const VectorXd &y_in, const VectorXd &z_in,
               const spline_endpoint_direction_3d &start_direction,
               const spline_endpoint_direction_3d &end_direction)
{
    const int n=x_in.size();
    if(n!=y_in.size() || n!=z_in.size())
        return;
    dimension=three;
    x=x_in;
    y=y_in;
    z=z_in;
    h=calculate_chord(x,y,z);
    VectorXd tmp0, tmp1, tmp2;
    if(!spline_arithmatic(h,x,start_direction.x_direction,end_direction.y_direction,tmp0)
            ||!spline_arithmatic(h,y,start_direction.x_direction,end_direction.y_direction,tmp1)
            ||!spline_arithmatic(h,z,start_direction.x_direction,end_direction.y_direction,tmp2))
        return;
    dif_2.resize(n,3);
    dif_2.col(0)=tmp0;
    dif_2.col(1)=tmp1;
    dif_2.col(2)=tmp2;
}

MatrixXd spline::get_dif_1()
{
    dif_1.resizeLike(dif_2);
    dif_1(0)=(y(1)-y(0))/h(0)+(2*dif_2(0)-dif_2(1))*h(0)/6;
    const int n=y.size();
    dif_1.block(1,0,n-1,1)=h.block(0,0,n-1,1).cwiseProduct(dif_2.block(0,0,n-1,1))/6
            +h.block(0,0,n-1,1).cwiseProduct(dif_2.block(1,0,n-1,1))/3
            +(y.block(1,0,n-1,1)-y.block(0,0,n-1,1)).cwiseQuotient(h.block(0,0,n-1,1));
    return(dif_1);
}

MatrixXd spline::get_coefficient()
{
    const int n=h.size();
    if(coefficient.rows()!=n)
        calculate_coefficient();
    return(coefficient);
}

VectorXd spline::get_length()
{
    const int n=h.size();
    if(coefficient.rows()!=n)
        calculate_coefficient();
    VectorXd tmp0, tmp1, tmp2, tmp3;
    for(int n1=0;n1<n;++n1)
    {
        tmp0.setLinSpaced(INTEGRATION_CONSTANT+1,0,h(n1));
        tmp1=coefficient(n1,2)*tmp0.cwiseProduct(tmp0).cwiseProduct(tmp0)
                +coefficient(n1,1)*tmp0.cwiseProduct(tmp0)
                +coefficient(n1,0)*tmp0;
        tmp2=coefficient(n1,5)*tmp0.cwiseProduct(tmp0).cwiseProduct(tmp0)
                +coefficient(n1,4)*tmp0.cwiseProduct(tmp0)
                +coefficient(n1,3)*tmp0;
        switch(dimension)
        {
        case spline::two:
            length(n1)=(tmp1.cwiseProduct(tmp1)+tmp2.cwiseProduct(tmp2)).cwiseSqrt().sum();
            break;
        case spline::three:
            tmp3=coefficient(n1,8)*tmp0.cwiseProduct(tmp0).cwiseProduct(tmp0)
                    +coefficient(n1,7)*tmp0.cwiseProduct(tmp0)
                    +coefficient(n1,6)*tmp0;
            length(n1)=(tmp1.cwiseProduct(tmp1)+tmp2.cwiseProduct(tmp2)+tmp3.cwiseProduct(tmp3)).cwiseSqrt().sum();
            break;
        }
    }
    return(length);
}

void spline::calculate_coefficient()
{
    const int n=h.size();
    MatrixXd tmp0, tmp1, tmp2;
    switch(dimension)
    {
    case spline::two:
        coefficient.resize(n,2*3);
        if(!spline_coefficient_arithmatic(h,x,dif_2.col(0),tmp0)
                ||!spline_coefficient_arithmatic(h,y,dif_2.col(1),tmp1))
            return;
        coefficient.leftCols<3>()=tmp0;
        coefficient.rightCols<3>()=tmp1;
        break;
    case spline::three:
        coefficient.resize(n,3*3);
        if(!spline_coefficient_arithmatic(h,x,dif_2.col(0),tmp0)
                ||!spline_coefficient_arithmatic(h,y,dif_2.col(1),tmp1)
                ||!spline_coefficient_arithmatic(h,z,dif_2.col(2),tmp2))
            return;
        coefficient.block(0,0,n,3)=tmp0;
        coefficient.block(0,3,n,3)=tmp1;
        coefficient.block(0,6,n,3)=tmp1;
        break;
    }
}

//core arithmatic
bool spline::spline_coefficient_arithmatic(const VectorXd &h, const VectorXd &y, const VectorXd &dif_2,
                                           MatrixXd &coeff)
{
    const int n=h.size();
    if(y.size()!=n+1
            ||dif_2.size()!=n+1)
        return(false);
    coeff.resize(n,3);
    coeff.col(0)=(dif_2.tail(n)-dif_2.head(n)).cwiseQuotient(6*h);
    coeff.col(1)=dif_2.head(n)/2;
    coeff.col(2)=(y.tail(n)-y.head(n)).cwiseQuotient(h)-(2*dif_2.head(n)+dif_2.tail(n)).cwiseProduct(h)/6;
    return(true);
}

bool spline::spline_arithmatic(const VectorXd &h, const VectorXd &y,
                               const spline_endpoint &start_point,
                               const spline_endpoint &end_point,
                               VectorXd &dif_2)
{
    const int n=y.size();
    if(n<3
            ||n-1!=h.size())
    {
        return(false);
    }
    VectorXd a(n-1),b(n-1),d(n);//for matrix construction
    dif_2.resize(n);
    switch(start_point.dif)
    {
        case spline_endpoint::first:
            a(0)=1;
            d(0)=6/h(0)*((y(1)-y(0))/h(0)-start_point.value);
            break;
        case spline_endpoint::second:
            a(0)=0;
            d(0)=2*start_point.value;
            break;
    }
    switch(end_point.dif)
    {
    case spline_endpoint::first:
        b(n-2)=1;
        d(n-1)=6/h(0)*(end_point.value-(y(n-1)-y(n-2))/h(n-2));
        break;
    case spline_endpoint::second:
        b(n-2)=0;
        d(n-1)=2*end_point.value;
        break;
    }
    d.block(1,0,n-2,1)=((y.block(2,0,n-2,1)-y.block(1,0,n-2,1)).cwiseQuotient(h.block(1,0,n-2,1))
                        -(y.block(1,0,n-2,1)-y.block(0,0,n-2,1)).cwiseQuotient(h.block(0,0,n-2,1))
                        ).cwiseQuotient(h.block(0,0,n-2,1)+h.block(1,0,n-2,1))*6;
    a.block(1,0,n-2,1)=h.block(1,0,n-2,1).cwiseQuotient(h.block(0,0,n-2,1)+h.block(1,0,n-2,1));
    b.block(0,0,n-2,1)=MatrixXd::Ones(n-2,1)-a.block(1,0,n-2,1);
    return(thomas_solver3(a, b, d, dif_2));
}

bool spline::thomas_solver3(const VectorXd &a, const VectorXd &b, const VectorXd &d, VectorXd &c)
{
    int n=d.size();
    if(n-1!=a.size() || n-1!=b.size())
        return(false);
    c.resize(n);
    VectorXd u(n-1),p(n-1),q(n-1);//for matrix solver
    u(0)=2;p(0)=0;q(0)=d(0);
    for(int i=1;i<n-1;++i)
    {
        if(0==u(i-1))
            return(false);
        p(i)=b(i)/u(i-1);
        u(i)=2-p(i)*a(i);
        q(i)=d(i)-p(i)*q(i-1);
    }
    c(n-1)=q(n-2)/u(n-2);c(0)=0;
    for(int i=n-2;i>0;--i)
    {
        c(i)=(q(i)-a(i)*c(i+1))/u(i);
    }
    return(true);
}

VectorXd spline::calculate_chord(const VectorXd &x, const VectorXd &y)
{
    const int n=x.size();
    VectorXd tmp1=x.tail(n-1)-x.head(n-1), tmp2=y.tail(n-1)-y.head(n-1);
    return((tmp1.cwiseProduct(tmp1)+tmp2.cwiseProduct(tmp2)).cwiseSqrt());
}

VectorXd spline::calculate_chord(const VectorXd &x, const VectorXd &y, const VectorXd &z)
{
    const int n=x.size();
    VectorXd tmp1=x.tail(n-1)-x.head(n-1),
            tmp2=y.tail(n-1)-y.head(n-1),
            tmp3=z.tail(n-1)-z.head(n-1);
    return((tmp1.cwiseProduct(tmp1)+tmp2.cwiseProduct(tmp2)+tmp3.cwiseProduct(tmp3)).cwiseSqrt());
}


spline_endpoint::spline_endpoint(order d, double v)
{
    dif=d;
    value=v;
}

spline_endpoint_direction_2d::spline_endpoint_direction_2d(const Vector2d &in):
    x_direction(spline_endpoint::first,in(0)),y_direction(spline_endpoint::first,in(1))
{
}

spline_endpoint_direction_2d::spline_endpoint_direction_2d(const spline_endpoint &x_dir,
                                                           const spline_endpoint &y_dir):
    x_direction(x_dir.dif,x_dir.value),y_direction(y_dir.dif,y_dir.value)
{
}

spline_endpoint_direction_3d::spline_endpoint_direction_3d(const Vector3d &in):
    x_direction(spline_endpoint::first,in(0)),
    y_direction(spline_endpoint::first,in(1)),
    z_direction(spline_endpoint::first,in(2))
{
}

spline_endpoint_direction_3d::spline_endpoint_direction_3d(const spline_endpoint &x_dir,
                                                           const spline_endpoint &y_dir,
                                                           const spline_endpoint &z_dir):
    x_direction(x_dir.dif,x_dir.value),
    y_direction(y_dir.dif,y_dir.value),
    z_direction(z_dir.dif,z_dir.value)
{
}

#endif //INTEGRATION_CONSTANT
