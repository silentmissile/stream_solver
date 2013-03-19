#ifndef INTEGRATION_CONSTANT
#define INTEGRATION_CONSTANT 20
#include "spline.h"

spline::spline()
{
}

//bool spline::spline3(const VectorXd &x, const VectorXd &y, VectorXd &dif_1) const
//{
//    const int n=x.rows();
//    if(n<3
//            ||n!=y.rows()
//            ||n!=dif_1.rows())
//    {
//        return(false);
//    }
//    VectorXd h=x.bottomRows(n-1)-x.topRows(n-1);
//    VectorXd d(n,1),c(n,1);
//    SparseMatrix<double> cal(n,n);

//    cal.reserve(VectorXd::Constant(n,3));
//    cal.insert(0,0)=2;
//    cal.insert(0,1)=0;
//    cal.insert(n-1,n-1)=2;
//    cal.insert(n-1,n-2)=0;
//    d(0)=0;
//    d(n-1)=0;
//    for(int i=1;i<n-1;++i)
//    {
//        cal.insert(i,i-1)=h(i-1)/(h(i-1)+h(i));
//        cal.insert(i,i+1)=h(i)/(h(i-1)+h(i));
//        cal.insert(i,i)=2;
//        d(i)=((y(i+1)-y(i))/h(i)-(y(i)-y(i-1))/h(i-1))*6/(h(i-1)+h(i));
//    }
//    cal.makeCompressed();

////    std::vector<Eigen::Triplet<double> > tripletList;
////    tripletList.push_back(Eigen::Triplet<double>(0,0,2));
////    tripletList.push_back(Eigen::Triplet<double>(0,1,0));
////    tripletList.push_back(Eigen::Triplet<double>(n-1,n-1,2));
////    tripletList.push_back(Eigen::Triplet<double>(n-1,n-2,0));
////    for(int i=1;i<n-1;++i)
////    {
////        tripletList.push_back(Eigen::Triplet<double>(i,i-1,h(i-1)/(h(i-1)+h(i))));
////        tripletList.push_back(Eigen::Triplet<double>(i,i+1,h(i)/(h(i-1)+h(i))));
////        tripletList.push_back(Eigen::Triplet<double>(i,i,2));
////        d(i)=((y(i+1)-y(i))/h(i)-(y(i)-y(i-1))/h(i-1))*6/(h(i-1)+h(i));
////    }
////    cal.setFromTriplets(tripletList.begin(), tripletList.end());
//    SimplicialLDLT<SparseMatrix<double> > solver;
//    solver.compute(cal);
//    if(solver.info()!=Success)
//    {
//      // decomposition failed
//        return(false);
//    }
//    c=solver.solve(d);
//    if(solver.info()!=Success)
//    {
//      // decomposition failed
//        return(false);
//    }
//    dif_1(0)=(y(1)-y(0))/h(0)-d(0)*h(0)/6;
//    //dif_1(n-1)=(y(n-1)-y(n-2))/h(n-2)-d(n-1)*h(n-2)/6;
////    for(int i=1;i<n;++i)
////    {
////        dif_1(i)=h(i-1)*c(i-1)/6+h(i-1)*c(i)/3+(y(i)-y(i-1))/h(i-1);
////    }
//    dif_1.block(1,0,n-1,1)=h.block(0,0,n-1,1).cwiseProduct(c.block(0,0,n-1,1))/6
//            +h.block(0,0,n-1,1).cwiseProduct(c.block(1,0,n-1,1))/3
//            +(y.block(1,0,n-1,1)-y.block(0,0,n-1,1)).cwiseQuotient(h.block(0,0,n-1,1));
//    return(true);
//}

bool spline::spline_parametric(const VectorXd &x, const VectorXd &y, VectorXd &dif_1)
{
    const int n=x.size();
    if(n<3 || n!=y.size())
        return(false);
    dif_1.resize(n);
    VectorXd h=x.tail(n-1)-x.head(n-1);
    if(spline_arithmatic(h,y,dif_1))
        return(true);
    else
        return(false);
}

bool spline::spline_parametric(const VectorXd &x, const VectorXd &y, MatrixXd &dif_1)
{
    const int n=x.size();
    if(n<3 || n!=y.size())
        return(false);
    dif_1.resize(n,2);
    VectorXd s(n-1);
    VectorXd tmpd0, tmpd1;
    tmpd0=x.tail(n-1)-x.head(n-1);
    tmpd1=y.tail(n-1)-y.head(n-1);
    s=(tmpd0.cwiseProduct(tmpd0)+tmpd1.cwiseProduct(tmpd1)).cwiseSqrt();
//    tmpd0.resize(n);
//    tmpd1.resize(n);
    if(!spline_arithmatic(s,x,tmpd0)
            ||!spline_arithmatic(s,y,tmpd1))
        return(false);
    dif_1.col(0)=tmpd0;
    dif_1.col(1)=tmpd1;
    return(true);
}

bool spline::spline_parametric(const VectorXd &x, const VectorXd &y,
                              MatrixXd &dif_1, MatrixXd &dif_2)
{
    const int n=x.size();
    if(n<3 || n!=y.size())
        return(false);
    dif_1.resize(n,2);
    dif_2.resize(n,2);
    VectorXd s(n-1);
    VectorXd tmpd0, tmpd1;
    tmpd0=x.tail(n-1)-x.head(n-1);
    tmpd1=y.tail(n-1)-y.head(n-1);
    s=(tmpd0.cwiseProduct(tmpd0)+tmpd1.cwiseProduct(tmpd1)).cwiseSqrt();
//    tmpd0.resize(n);
//    tmpd1.resize(n);
    VectorXd tmpd2(n), tmpd3(n);
    if(!spline_arithmatic(s,x,tmpd0,tmpd1)
            ||!spline_arithmatic(s,y,tmpd2,tmpd3))
        return(false);
    dif_1.col(0)=tmpd0;
    dif_1.col(1)=tmpd1;
    dif_2.col(0)=tmpd2;
    dif_2.col(1)=tmpd3;
    return(true);
}

bool spline::spline_parametric(const VectorXd &x, const VectorXd &y,
                               const Vector2d &start_direction, const Vector2d &end_direction,
                               MatrixXd &dif_1)
{
    const int n=x.size();
    if(n<3 || n!=y.size())
        return(false);
    dif_1.resize(n,2);
    VectorXd s(n-1);
    VectorXd tmpd0, tmpd1;
    tmpd0=x.tail(n-1)-x.head(n-1);
    tmpd1=y.tail(n-1)-y.head(n-1);
    s=(tmpd0.cwiseProduct(tmpd0)+tmpd1.cwiseProduct(tmpd1)).cwiseSqrt();
//    tmpd0.resize(n);
//    tmpd1.resize(n);
    if(!spline_arithmatic(s,x,start_direction(0),end_direction(0),tmpd0)
            ||!spline_arithmatic(s,y,start_direction(1),end_direction(1),tmpd1))
        return(false);
    dif_1.col(0)=tmpd0;
    dif_1.col(1)=tmpd1;
    return(true);
}

bool spline::spline_parametric(const VectorXd &x, const VectorXd &y,
                               const Vector2d &start_direction, const Vector2d &end_direction,
                               MatrixXd &dif_1, MatrixXd &dif_2)
{
    const int n=x.size();
    if(n<3 || n!=y.size())
        return(false);
    dif_1.resize(n,2);
    dif_2.resize(n,2);
    VectorXd s(n-1);
    VectorXd tmpd0, tmpd1;
    tmpd0=x.tail(n-1)-x.head(n-1);
    tmpd1=y.tail(n-1)-y.head(n-1);
    s=(tmpd0.cwiseProduct(tmpd0)+tmpd1.cwiseProduct(tmpd1)).cwiseSqrt();
//    tmpd0.resize(n);
//    tmpd1.resize(n);
    VectorXd tmpd2(n), tmpd3(n);
    if(!spline_arithmatic(s,x,start_direction(0),end_direction(0),tmpd0,tmpd2)
            ||!spline_arithmatic(s,y,start_direction(1),end_direction(1),tmpd1,tmpd3))
        return(false);
    dif_1.col(0)=tmpd0;
    dif_1.col(1)=tmpd1;
    dif_2.col(0)=tmpd2;
    dif_2.col(1)=tmpd3;
    return(true);
}

bool spline::spline_parametric(const VectorXd &x, const VectorXd &y,
                       const Vector2d &start_direction, const Vector2d &end_direction,
                       MatrixXd &dif_1, MatrixXd &dif_2, VectorXd &len)
{
    const int n=x.size();
    if(n<3 || n!=y.size())
        return(false);
    dif_1.resize(n,2);
    dif_2.resize(n,2);
    VectorXd s(n-1);
    VectorXd tmpd0, tmpd1;
    tmpd0=x.tail(n-1)-x.head(n-1);
    tmpd1=y.tail(n-1)-y.head(n-1);
    s=(tmpd0.cwiseProduct(tmpd0)+tmpd1.cwiseProduct(tmpd1)).cwiseSqrt();
//    tmpd0.resize(n);
//    tmpd1.resize(n);
    VectorXd tmpd2(n), tmpd3(n);
    if(!spline_arithmatic(s,x,start_direction(0),end_direction(0),tmpd0,tmpd2)
            ||!spline_arithmatic(s,y,start_direction(1),end_direction(1),tmpd1,tmpd3))
        return(false);
    dif_1.col(0)=tmpd0;
    dif_1.col(1)=tmpd1;
    dif_2.col(0)=tmpd2;
    dif_2.col(1)=tmpd3;
    MatrixXd coef(n-1,8);
//    VectorXd s0(n);
//    s0(0)=0;
//    for(int k=1;k<n;++k)
//    {
//        s0(k)=s0(k-1)+s(k-1);
//    }
    coef.col(0)=x.head(n-1).cwiseQuotient(s);
    coef.col(1)=(1/3*dif_2.col(0).head(n-1)-1/6*dif_2.col(0).tail(n-1)).cwiseProduct(s)
            +(x.tail(n-1)-x.head(n-1)).cwiseQuotient(s);
    coef.col(2)=0.5*dif_2.col(0).cwiseProduct(s);
    coef.col(3)=1/6*(dif_2.col(0).tail(n-1)-dif_2.col(0).head(n-1));
    coef.col(4)=y.head(n-1).cwiseQuotient(s);
    coef.col(5)=(1/3*dif_2.col(0).head(n-1)-1/6*dif_2.col(0).tail(n-1)).cwiseProduct(s)
            +(y.tail(n-1)-y.head(n-1)).cwiseQuotient(s);
    coef.col(6)=coef.col(2);
    coef.col(7)=coef.col(3);
    len.resizeLike(s);
//    tmpd0.resize(INTEGRATION_CONSTANT+1);
//    tmpd1.resize(INTEGRATION_CONSTANT+1);
//    tmpd2.resize(INTEGRATION_CONSTANT+1);
    for(int n1=0;n1<n-1;++n1)
    {
        tmpd0.setLinSpaced(INTEGRATION_CONSTANT+1,0,s(n1));
        tmpd1=coef(n1,3)*tmpd0.cwiseProduct(tmpd0).cwiseProduct(tmpd0)
                +coef(n1,2)*tmpd0.cwiseProduct(tmpd0)
                +coef(n1,1)*tmpd0;
        tmpd2=coef(n1,7)*tmpd0.cwiseProduct(tmpd0).cwiseProduct(tmpd0)
                +coef(n1,6)*tmpd0.cwiseProduct(tmpd0)
                +coef(n1,5)*tmpd0;
        len(n1)=(tmpd1.cwiseProduct(tmpd1)+tmpd2.cwiseProduct(tmpd2)).cwiseSqrt().sum();
    }
    return(true);
}

bool spline::spline_parametric(const VectorXd &x, const VectorXd &y, const VectorXd &z,
                               MatrixXd &dif_1)
{
    const int n=x.size();
    if(n<3
            ||n!=y.size()
            ||n!=z.size())
        return(false);
    dif_1.resize(n,3);
    VectorXd s(n-1);
    VectorXd tmpd0, tmpd1, tmpd2;
    tmpd0=x.tail(n-1)-x.head(n-1);
    tmpd1=y.tail(n-1)-y.head(n-1);
    tmpd2=z.tail(n-1)-z.head(n-1);
    s=(tmpd0.cwiseProduct(tmpd0)+tmpd1.cwiseProduct(tmpd1)+tmpd2.cwiseProduct(tmpd2)).cwiseSqrt();
//    tmpd0.resize(n);
//    tmpd1.resize(n);
//    tmpd2.resize(n);
    if(!spline_arithmatic(s,x,tmpd0)
            ||!spline_arithmatic(s,y,tmpd1)
            ||!spline_arithmatic(s,z,tmpd2))
        return(false);
    dif_1.col(0)=tmpd0;
    dif_1.col(1)=tmpd1;
    dif_1.col(2)=tmpd2;
    return(true);
}

bool spline::spline_parametric(const VectorXd &x, const VectorXd &y, const VectorXd &z,
                               MatrixXd &dif_1, MatrixXd &dif_2)
{
    const int n=x.size();
    if(n<3
            ||n!=y.size()
            ||n!=z.size())
        return(false);
    dif_1.resize(n,3);
    VectorXd s(n-1);
    VectorXd tmpd0, tmpd1, tmpd2;
    tmpd0=x.tail(n-1)-x.head(n-1);
    tmpd1=y.tail(n-1)-y.head(n-1);
    tmpd2=z.tail(n-1)-z.head(n-1);
    s=(tmpd0.cwiseProduct(tmpd0)+tmpd1.cwiseProduct(tmpd1)+tmpd2.cwiseProduct(tmpd2)).cwiseSqrt();
//    tmpd0.resize(n);
//    tmpd1.resize(n);
//    tmpd2.resize(n);
    VectorXd tmpd3(n), tmpd4(n), tmpd5(n);
    if(!spline_arithmatic(s,x,tmpd0,tmpd3)
            ||!spline_arithmatic(s,y,tmpd1,tmpd4)
            ||!spline_arithmatic(s,z,tmpd2,tmpd5))
        return(false);
    dif_1.col(0)=tmpd0;
    dif_1.col(1)=tmpd1;
    dif_1.col(2)=tmpd2;
    dif_2.col(0)=tmpd3;
    dif_2.col(1)=tmpd4;
    dif_2.col(2)=tmpd5;
    return(true);
}

bool spline::spline_parametric(const VectorXd &x, const VectorXd &y, const VectorXd &z,
                               const Vector3d &start_direction, const Vector3d &end_direction,
                               MatrixXd &dif_1)
{
    const int n=x.size();
    if(n<3
            ||n!=y.size()
            ||n!=z.size())
        return(false);
    dif_1.resize(n,3);
    VectorXd s(n-1);
    VectorXd tmpd0, tmpd1, tmpd2;
    tmpd0=x.tail(n-1)-x.head(n-1);
    tmpd1=y.tail(n-1)-y.head(n-1);
    tmpd2=z.tail(n-1)-z.head(n-1);
    s=(tmpd0.cwiseProduct(tmpd0)+tmpd1.cwiseProduct(tmpd1)+tmpd2.cwiseProduct(tmpd2)).cwiseSqrt();
//    tmpd0.resize(n);
//    tmpd1.resize(n);
//    tmpd2.resize(n);
    if(!spline_arithmatic(s,x,start_direction(0),end_direction(0),tmpd0)
            ||!spline_arithmatic(s,y,start_direction(1),end_direction(1),tmpd1)
            ||!spline_arithmatic(s,z,start_direction(2),end_direction(2),tmpd2))
        return(false);
    dif_1.col(0)=tmpd0;
    dif_1.col(1)=tmpd1;
    dif_1.col(2)=tmpd2;
    return(true);
}

bool spline::spline_parametric(const VectorXd &x, const VectorXd &y, const VectorXd &z,
                                  const Vector3d &start_direction, const Vector3d &end_direction,
                                  MatrixXd &dif_1, MatrixXd &dif_2)
{
    const int n=x.size();
    if(n<3
            ||n!=y.size()
            ||n!=z.size())
        return(false);
    dif_1.resize(n,3);
    VectorXd s(n-1);
    VectorXd tmpd0, tmpd1, tmpd2;
    tmpd0=x.tail(n-1)-x.head(n-1);
    tmpd1=y.tail(n-1)-y.head(n-1);
    tmpd2=z.tail(n-1)-z.head(n-1);
    s=(tmpd0.cwiseProduct(tmpd0)+tmpd1.cwiseProduct(tmpd1)+tmpd2.cwiseProduct(tmpd2)).cwiseSqrt();
//    tmpd0.resize(n);
//    tmpd1.resize(n);
//    tmpd2.resize(n);
    VectorXd tmpd3(n), tmpd4(n), tmpd5(n);
    if(!spline_arithmatic(s,x,start_direction(0),end_direction(0),tmpd0,tmpd3)
            ||!spline_arithmatic(s,y,start_direction(1),end_direction(1),tmpd1,tmpd4)
            ||!spline_arithmatic(s,z,start_direction(2),end_direction(2),tmpd2,tmpd5))
        return(false);
    dif_1.col(0)=tmpd0;
    dif_1.col(1)=tmpd1;
    dif_1.col(2)=tmpd2;
    dif_2.col(0)=tmpd3;
    dif_2.col(1)=tmpd4;
    dif_2.col(2)=tmpd5;
    return(true);
}

bool spline::spline_parametric(const VectorXd &x, const VectorXd &y, const VectorXd &z,
                                  const Vector3d &start_direction, const Vector3d &end_direction,
                                  MatrixXd &dif_1, MatrixXd &dif_2, VectorXd &len)
{
    const int n=x.size();
    if(n<3
            ||n!=y.size()
            ||n!=z.size())
        return(false);
    dif_1.resize(n,3);
    VectorXd s(n-1);
    VectorXd tmpd0, tmpd1, tmpd2;
    tmpd0=x.tail(n-1)-x.head(n-1);
    tmpd1=y.tail(n-1)-y.head(n-1);
    tmpd2=z.tail(n-1)-z.head(n-1);
    s=(tmpd0.cwiseProduct(tmpd0)+tmpd1.cwiseProduct(tmpd1)+tmpd2.cwiseProduct(tmpd2)).cwiseSqrt();
//    tmpd0.resize(n);
//    tmpd1.resize(n);
//    tmpd2.resize(n);
    VectorXd tmpd3(n), tmpd4(n), tmpd5(n);
    if(!spline_arithmatic(s,x,start_direction(0),end_direction(0),tmpd0,tmpd3)
            ||!spline_arithmatic(s,y,start_direction(1),end_direction(1),tmpd1,tmpd4)
            ||!spline_arithmatic(s,z,start_direction(2),end_direction(2),tmpd2,tmpd5))
        return(false);
    dif_1.col(0)=tmpd0;
    dif_1.col(1)=tmpd1;
    dif_1.col(2)=tmpd2;
    dif_2.col(0)=tmpd3;
    dif_2.col(1)=tmpd4;
    dif_2.col(2)=tmpd5;
    MatrixXd coef(n-1,12);
//    VectorXd s0(n);
//    s0(0)=0;
//    for(int k=1;k<n;++k)
//    {
//        s0(k)=s0(k-1)+s(k-1);
//    }
    coef.col(0)=x.head(n-1).cwiseQuotient(s);
    coef.col(1)=(1/3*dif_2.col(0).head(n-1)-1/6*dif_2.col(0).tail(n-1)).cwiseProduct(s)
            +(x.tail(n-1)-x.head(n-1)).cwiseQuotient(s);
    coef.col(2)=0.5*dif_2.col(0).cwiseProduct(s);
    coef.col(3)=1/6*(dif_2.col(0).tail(n-1)-dif_2.col(0).head(n-1));
    coef.col(4)=y.head(n-1).cwiseQuotient(s);
    coef.col(5)=(1/3*dif_2.col(0).head(n-1)-1/6*dif_2.col(0).tail(n-1)).cwiseProduct(s)
            +(y.tail(n-1)-y.head(n-1)).cwiseQuotient(s);
    coef.col(6)=coef.col(2);
    coef.col(7)=coef.col(3);
    coef.col(8)=z.head(n-1).cwiseQuotient(s);
    coef.col(9)=(1/3*dif_2.col(0).head(n-1)-1/6*dif_2.col(0).tail(n-1)).cwiseProduct(s)
            +(z.tail(n-1)-z.head(n-1)).cwiseQuotient(s);
    coef.col(10)=coef.col(2);
    coef.col(11)=coef.col(3);
    len.resizeLike(s);
//    tmpd0.resize(INTEGRATION_CONSTANT+1);
//    tmpd1.resize(INTEGRATION_CONSTANT+1);
//    tmpd2.resize(INTEGRATION_CONSTANT+1);
//    tmpd3.resize(INTEGRATION_CONSTANT+1);
    for(int n1=0;n1<n-1;++n1)
    {
        tmpd0.setLinSpaced(INTEGRATION_CONSTANT+1,0,s(n1));
        tmpd1=coef(n1,3)*tmpd0.cwiseProduct(tmpd0).cwiseProduct(tmpd0)
                +coef(n1,2)*tmpd0.cwiseProduct(tmpd0)
                +coef(n1,1)*tmpd0;
        tmpd2=coef(n1,7)*tmpd0.cwiseProduct(tmpd0).cwiseProduct(tmpd0)
                +coef(n1,6)*tmpd0.cwiseProduct(tmpd0)
                +coef(n1,5)*tmpd0;
        tmpd3=coef(n1,11)*tmpd0.cwiseProduct(tmpd0).cwiseProduct(tmpd0)
                +coef(n1,10)*tmpd0.cwiseProduct(tmpd0)
                +coef(n1,9)*tmpd0;
        len(n1)=(tmpd1.cwiseProduct(tmpd1)+tmpd2.cwiseProduct(tmpd2)+tmpd3.cwiseProduct(tmpd3)).cwiseSqrt().sum();
    }
    return(true);
}

//core arithmatic
bool spline::spline_arithmatic(const VectorXd &h, const VectorXd &y, VectorXd &dif_1)
{
    const int n=y.size();
    if(n<3
            ||n-1!=h.size())
    {
        return(false);
    }
    dif_1.resize(n);
    VectorXd a(n-1),b(n-1),d(n),c(n);//for matrix construction
    //vector c is 2nd order differential
    d(0)=0;d(n-1)=0;
    d.block(1,0,n-2,1)=((y.block(2,0,n-2,1)-y.block(1,0,n-2,1)).cwiseQuotient(h.block(1,0,n-2,1))
                        -(y.block(1,0,n-2,1)-y.block(0,0,n-2,1)).cwiseQuotient(h.block(0,0,n-2,1))
                        ).cwiseQuotient(h.block(0,0,n-2,1)+h.block(1,0,n-2,1))*6;
    a(0)=0;b(n-2)=0;
    a.block(1,0,n-2,1)=h.block(1,0,n-2,1).cwiseQuotient(h.block(0,0,n-2,1)+h.block(1,0,n-2,1));
    b.block(0,0,n-2,1)=MatrixXd::Ones(n-2,1)-a.block(1,0,n-2,1);
    if(!thomas_solver3(a, b, d, c))
        return(false);
    dif_1(0)=(y(1)-y(0))/h(0)-d(0)*h(0)/6;
    dif_1.block(1,0,n-1,1)=h.block(0,0,n-1,1).cwiseProduct(c.block(0,0,n-1,1))/6
            +h.block(0,0,n-1,1).cwiseProduct(c.block(1,0,n-1,1))/3
            +(y.block(1,0,n-1,1)-y.block(0,0,n-1,1)).cwiseQuotient(h.block(0,0,n-1,1));
    return(true);
}

bool spline::spline_arithmatic(const VectorXd &h, const VectorXd &y,
                               const double &start_dif_1, const double &end_dif_1,
                               VectorXd &dif_1)
{
    const int n=y.size();
    if(n<3
            ||n-1!=h.size())
    {
        return(false);
    }
    dif_1.resize(n);
    VectorXd a(n-1),b(n-1),d(n),c(n);//for matrix construction
    //vector c is 2nd order differential
    d(0)=6/h(0)*((y(1)-y(0))/h(0)-start_dif_1);
    d(n-1)=6/h(0)*(end_dif_1-(y(n-1)-y(n-2))/h(n-2));
    d.block(1,0,n-2,1)=((y.block(2,0,n-2,1)-y.block(1,0,n-2,1)).cwiseQuotient(h.block(1,0,n-2,1))
                        -(y.block(1,0,n-2,1)-y.block(0,0,n-2,1)).cwiseQuotient(h.block(0,0,n-2,1))
                        ).cwiseQuotient(h.block(0,0,n-2,1)+h.block(1,0,n-2,1))*6;
    a(0)=1;
    b(n-2)=1;
    a.block(1,0,n-2,1)=h.block(1,0,n-2,1).cwiseQuotient(h.block(0,0,n-2,1)+h.block(1,0,n-2,1));
    b.block(0,0,n-2,1)=MatrixXd::Ones(n-2,1)-a.block(1,0,n-2,1);
    if(!thomas_solver3(a, b, d, c))
        return(false);
    dif_1(0)=(y(1)-y(0))/h(0)-d(0)*h(0)/6;
    dif_1.block(1,0,n-1,1)=h.block(0,0,n-1,1).cwiseProduct(c.block(0,0,n-1,1))/6
            +h.block(0,0,n-1,1).cwiseProduct(c.block(1,0,n-1,1))/3
            +(y.block(1,0,n-1,1)-y.block(0,0,n-1,1)).cwiseQuotient(h.block(0,0,n-1,1));
    return(true);
}

bool spline::spline_arithmatic(const VectorXd &h, const VectorXd &y, VectorXd &dif_1, VectorXd &dif_2)
{
    const int n=y.size();
    if(n<3
            ||n-1!=h.size())
    {
        return(false);
    }
    dif_1.resize(n);
    dif_2.resize(n);
    VectorXd a(n-1),b(n-1),d(n);//for matrix construction
    d(0)=0;d(n-1)=0;
    d.block(1,0,n-2,1)=((y.block(2,0,n-2,1)-y.block(1,0,n-2,1)).cwiseQuotient(h.block(1,0,n-2,1))
                        -(y.block(1,0,n-2,1)-y.block(0,0,n-2,1)).cwiseQuotient(h.block(0,0,n-2,1))
                        ).cwiseQuotient(h.block(0,0,n-2,1)+h.block(1,0,n-2,1))*6;
    a(0)=0;b(n-2)=0;
    a.block(1,0,n-2,1)=h.block(1,0,n-2,1).cwiseQuotient(h.block(0,0,n-2,1)+h.block(1,0,n-2,1));
    b.block(0,0,n-2,1)=MatrixXd::Ones(n-2,1)-a.block(1,0,n-2,1);
    if(!thomas_solver3(a, b, d, dif_2))
        return(false);
    dif_1(0)=(y(1)-y(0))/h(0)-d(0)*h(0)/6;
    dif_1.block(1,0,n-1,1)=h.block(0,0,n-1,1).cwiseProduct(dif_2.block(0,0,n-1,1))/6
            +h.block(0,0,n-1,1).cwiseProduct(dif_2.block(1,0,n-1,1))/3
            +(y.block(1,0,n-1,1)-y.block(0,0,n-1,1)).cwiseQuotient(h.block(0,0,n-1,1));
    return(true);
}

bool spline::spline_arithmatic(const VectorXd &h, const VectorXd &y,
                               const double &start_dif_1, const double &end_dif_1,
                               VectorXd &dif_1, VectorXd &dif_2)
{
    const int n=y.size();
    if(n<3
            ||n-1!=h.size())
    {
        return(false);
    }
    VectorXd a(n-1),b(n-1),d(n);//for matrix construction
    dif_1.resize(n);
    dif_2.resize(n);
    d(0)=6/h(0)*((y(1)-y(0))/h(0)-start_dif_1);
    d(n-1)=6/h(0)*(end_dif_1-(y(n-1)-y(n-2))/h(n-2));
    d.block(1,0,n-2,1)=((y.block(2,0,n-2,1)-y.block(1,0,n-2,1)).cwiseQuotient(h.block(1,0,n-2,1))
                        -(y.block(1,0,n-2,1)-y.block(0,0,n-2,1)).cwiseQuotient(h.block(0,0,n-2,1))
                        ).cwiseQuotient(h.block(0,0,n-2,1)+h.block(1,0,n-2,1))*6;
    a(0)=1;
    b(n-2)=1;
    a.block(1,0,n-2,1)=h.block(1,0,n-2,1).cwiseQuotient(h.block(0,0,n-2,1)+h.block(1,0,n-2,1));
    b.block(0,0,n-2,1)=MatrixXd::Ones(n-2,1)-a.block(1,0,n-2,1);
    if(!thomas_solver3(a, b, d, dif_2))
        return(false);
    dif_1(0)=(y(1)-y(0))/h(0)-d(0)*h(0)/6;
    dif_1.block(1,0,n-1,1)=h.block(0,0,n-1,1).cwiseProduct(dif_2.block(0,0,n-1,1))/6
            +h.block(0,0,n-1,1).cwiseProduct(dif_2.block(1,0,n-1,1))/3
            +(y.block(1,0,n-1,1)-y.block(0,0,n-1,1)).cwiseQuotient(h.block(0,0,n-1,1));
    return(true);
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

#endif //INTEGRATION_CONSTANT
