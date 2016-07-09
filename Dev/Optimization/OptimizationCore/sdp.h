#ifndef SDP_H
#define SDP_H
#include "optimizationcore_global.h"
#include <armadillo>
#include <string>
namespace  Optimization {
class OPTIMIZATIONCORESHARED_EXPORT SDP
{
public:
    SDP();
    virtual ~SDP();
    void setC(const std::vector<arma::mat>&);
    void setAs(const std::vector<std::vector<arma::mat>>&);
    void setb(const arma::vec&);
    bool init(void);
    bool solve(void);
    void debug_prob(char*);
    void gety(arma::vec &y);
    void getX(arma::sp_mat& X);
    void getZ(arma::sp_mat& Z);
    virtual std::string info()const;
protected:

private:
    int n_,k_;
    double* b_;
    double *y_;
    double pobj_,dobj_;
    int ret_;
//pointer to private type
    void* C_;
    void* constraints_;
    void* X_;
    void* Z_;
};
}

#endif // SDP_H
