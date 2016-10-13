#include "jrcsaopt.h"
namespace JRCS {
JRCSAOPT::JRCSAOPT():JRCSAONI()
{

}
//constrain the alpha
//according to the input labels
//so that each patch in previous over-segmentation will not be split
void JRCSAOPT::alpha_operation(int i)
{
    arma::uvec& input_label = *vll_ptrlst_[i];
    arma::fmat& alpha = *alpha_ptrlst_[i];
    arma::uword lmax = input_label.max();
    arma::uword lmin = input_label.min();
    arma::fvec global_obj_prob_mean(obj_num_);
    std::vector<arma::uvec> oidxs(obj_num_);
    arma::frowvec global_alpha_mean = arma::mean(alpha);
    assert(global_alpha_mean.is_finite());
    #pragma omp parallel for
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        oidxs[o] = arma::find( obj_label_ == (o+1) );
        global_obj_prob_mean(o) = arma::accu( global_alpha_mean(oidxs[o]) );
    }
//    std::cerr<<"global_obj_prob_mean:"<<global_obj_prob_mean<<std::endl;
    while( lmin <= lmax )
    {
        arma::uvec current_patch_index = arma::find( lmin == input_label );
        ++lmin;
        if(current_patch_index.empty())continue;
        arma::fmat current_alpha = alpha.rows(current_patch_index);
        arma::frowvec current_alpha_mean = arma::mean( current_alpha );
        arma::fvec current_obj_prob_mean( obj_num_ );
       // #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            current_obj_prob_mean(o) = arma::accu( current_alpha_mean(oidxs[o]) );
            current_obj_prob_mean(o) /= global_obj_prob_mean(o);
        }
        current_obj_prob_mean += beta_; //add a small number to prevent underflow
        double sum = arma::accu(current_obj_prob_mean);
        assert(sum > 0);
        if( 0 != sum )current_obj_prob_mean /= sum;
        assert( current_obj_prob_mean.is_finite() );
//        std::cerr<<"d F"<<i<<"L"<<lmin<<" mean_obj_prob:"<<current_obj_prob_mean.t()<<std::endl;
        #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            alpha(current_patch_index,oidxs[o]) *= current_obj_prob_mean(o);
        }
    }
    //normalise alpha
    arma::fvec alpha_rowsum = arma::sum(alpha,1);
    alpha.each_col() /= alpha_rowsum;
}
}
