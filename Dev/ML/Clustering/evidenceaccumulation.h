#ifndef EVIDENCEACCUMULATION_H
#define EVIDENCEACCUMULATION_H
#include "clustering_global.h"
#include <armadillo>
#include <queue>
#include <functional>
#include "configure.h"
namespace Clustering{
class CLUSTERINGSHARED_EXPORT PEAC
{
public:
    typedef struct{
        double alpha_;
        double beta_;
        arma::uword gamma_;
    }DY;
    typedef struct Triplet{
        arma::uword index_;
        DY value_;
        double prior_;
        friend bool operator < (const struct Triplet &t1,const struct Triplet &t2){
                return t1.prior_ < t2.prior_;
            }
        friend bool operator > (const struct Triplet &t1,const struct Triplet &t2){
                return t1.prior_ > t2.prior_;
            }
    }Triplet;
    typedef struct {
        arma::uword index_;
        std::vector<arma::uword> indices_;
    }PNode;
    typedef struct {
        arma::uword index_;
        std::vector<double> value_;
    }ANode;
    PEAC();
    bool configure(Config::Ptr);
    void compute(
            const arma::umat& E,
            const arma::umat& P,
            arma::uvec& y
            );
protected:
    void compute();
    void computeCandN();
    void initPList();
    void initY();
    void initA();
    void initG();
    void initPQ();
    void computeD(); // Compute Update Direction
    void computeStep(); // Compute Step Size
    void computeY();    // Update Y
    void computeA();    // Update A
    void computeG();    // Update G
    void updatePrior(); // Update Prior
    void projectY(arma::uvec&);    // project Y
private:
    arma::umat iE_;
    arma::umat iP_;
    std::vector<Triplet> pq_;
    DY bestDY_;
    std::vector<PNode> PList_;
    std::shared_ptr<std::vector<ANode>> oldA_;
    std::shared_ptr<std::vector<ANode>> newA_;
    arma::sp_mat C_;
    double N_;
    arma::mat A_;
    arma::mat G_;
    arma::mat Y_;
    double tol_;
    arma::uword max_iter_;
};
}
#endif // EVIDENCEACCUMULATION_H
