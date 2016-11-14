#ifndef HKS_H
#define HKS_H
#include "featurecore_global.h"
#include "common.h"
namespace Feature{
template<typename Mesh>
class HKS
{
public:
    typedef std::shared_ptr<arma::mat> MatPtr;
    typedef std::vector<MatPtr> MatPtrLst;
    typedef std::vector<arma::uvec> LabelLst;
    HKS();
    bool configure(Config::Ptr);
    void extract(const typename MeshBundle<Mesh>::Ptr,arma::mat&);
    void extract(const typename MeshBundle<Mesh>::PtrList,MatPtrLst&);
    void extract(const typename MeshBundle<Mesh>::Ptr,const arma::uvec&,arma::mat&);
    void extract(const typename MeshBundle<Mesh>::PtrList,const LabelLst&,MatPtrLst&);
protected:
    double d_scale_;
    double convex_scale_;
    inline double distanceAffinity(double x1, double y1, double x2, double y2,double scale)
    {
        return -((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))/scale;
    }
    template<typename vec>
    inline double vecAffinity(const vec& x,const vec& y,double scale)
    {
        vec dif = x-y;
        return - arma::dot(dif,dif)/scale;
    }
    template<typename vec>
    inline double convexity(const vec& p0,const vec&n0,const vec&p1,const vec& n1,double eps)
    {
        eps = eps > 0?100.0*eps:std::numeric_limits<double>::epsilon();
        vec dir = arma::normalise(p1-p0);
        double dif = std::acos(arma::dot(n0,dir)) - std::acos(arma::dot(n1,dir));
        if(dif<-eps)return 0.0;
        else if(dif>eps)return 1.0;
        else return 0.5*dif/eps+0.5;
    }
    void constructL(const typename MeshBundle<Mesh>::Ptr);
    void constructL(const typename MeshBundle<Mesh>::Ptr,const arma::uvec&);
    void decomposeL();
    void computeHKS(arma::mat&);
private:
    arma::sp_mat W_;
    arma::vec lambda_;
    arma::mat eig_vec_;
};
}
#endif // HKS_H
