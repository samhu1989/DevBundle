#ifndef EVIDENCEACCUMULATION_H
#define EVIDENCEACCUMULATION_H
#include "clustering_global.h"
#include <armadillo>
#include <queue>
#include <functional>
#include "configure.h"
namespace Clustering{
//Consensus Clustering Using Partial Evidence Accumulation
class CLUSTERINGSHARED_EXPORT PEAC
{
public:
    typedef struct{
        arma::uword alpha_;
        arma::uword beta_;
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
    void initY_Random();
    void initY_RandIndex();
    void initA();
    void initG();
    void initPQ();
    void getBestD();
    void computeStep(); // Compute Step Size
    void computeY();    // Update Y
    void computeA();    // Update A
    void computeG();    // Update G
    void computePrior(Triplet&);
    void updatePrior(); // Update Prior
    double computeObj();
    void projectY(arma::uvec&);    // project Y
private:
    arma::umat iE_;
    arma::umat iP_;
    std::vector<Triplet> pq_;
    DY bestDY_;
    arma::uvec Pgamma_;
    std::shared_ptr<arma::sp_mat> oldA_;
    std::shared_ptr<arma::sp_mat> newA_;
    std::shared_ptr<arma::mat> newY_;
    std::shared_ptr<arma::mat> oldY_;
    arma::sp_mat C_;
    double N_;
    arma::mat G_;
    double step_;
    double tol_;
    arma::uword k_;
    arma::uword max_iter_;
    arma::uword init_;
    std::mt19937 rand_engine_;  // Mersenne twister random number engine
    double last_obj_;
    std::uniform_real_distribution<double> distr_;
};
}
#endif // EVIDENCEACCUMULATION_H
