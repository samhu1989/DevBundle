#include "jrcsaopt.h"
namespace JRCS {
JRCSAOPT::JRCSAOPT():JRCSAONI()
{

}
//constrain the alpha
//according to the input labels
//so that each patch in previous over-segmentation will not be split
//void JRCSAOPT::alpha_operation(int i)
//{
//    arma::uvec& input_label = *vll_ptrlst_[i];
//    arma::mat& alpha = *alpha_ptrlst_[i];
//    arma::uword lmax = input_label.max();
//    arma::uword lmin = input_label.min();
//    arma::fvec global_obj_prob_mean(obj_num_);
//    std::vector<arma::uvec> oidxs(obj_num_);
//    arma::rowvec global_alpha_mean = arma::mean(alpha);
//    assert(global_alpha_mean.is_finite());
//    #pragma omp parallel for
//    for(int o = 0 ; o < obj_num_ ; ++o )
//    {
//        oidxs[o] = arma::find( obj_label_ == (o+1) );
//        global_obj_prob_mean(o) = arma::accu( global_alpha_mean(oidxs[o]) );
//    }
////    std::cerr<<"global_obj_prob_mean:"<<global_obj_prob_mean<<std::endl;
//    while( lmin <= lmax )
//    {
//        arma::uvec current_patch_index = arma::find( lmin == input_label );
//        ++lmin;
//        if(current_patch_index.empty())continue;
//        arma::mat current_alpha = alpha.rows(current_patch_index);
//        arma::rowvec current_alpha_mean = arma::mean( current_alpha );
//        arma::fvec current_obj_prob_mean( obj_num_ );
//        #pragma omp parallel for
//        for(int o = 0 ; o < obj_num_ ; ++o )
//        {
//            current_obj_prob_mean(o) = arma::accu( current_alpha_mean(oidxs[o]) );
//            if( 0 != global_obj_prob_mean(o) )current_obj_prob_mean(o) /= global_obj_prob_mean(o);
//        }
//        double sum = arma::accu(current_obj_prob_mean);
////        std::cerr<<"sum:"<<sum<<std::endl;
//        assert(std::isfinite(sum));
//        if( 0 != sum )current_obj_prob_mean /= sum;
//        assert( current_obj_prob_mean.is_finite() );
////        std::cerr<<"d F"<<i<<"L"<<lmin<<" mean_obj_prob:"<<current_obj_prob_mean.t()<<std::endl;
//        #pragma omp parallel for
//        for(int o = 0 ; o < obj_num_ ; ++o )
//        {
//            alpha(current_patch_index,oidxs[o]) *= current_obj_prob_mean(o);
//        }
//    }
//    //normalise alpha
//    arma::vec alpha_rowsum = arma::sum(alpha,1);
//    alpha.each_col() /= alpha_rowsum;
//}
void JRCSAOPT::prepare_alpha_operation(int i)
{
//    std::cerr<<"alpha prepare"<<std::endl;
//    assert(false);
    if(iter_count_<=max_init_iter_)
    {
        if( i == init_alpha_value.size() )
        {
            arma::mat& alpha = *alpha_ptrlst_[i];
            init_alpha_value.emplace_back(
                        alpha.memptr(),
                        alpha.n_rows,
                        alpha.n_cols,
                        true,
                        true
                        );
        }
    }
    else if( iter_count_ == (max_init_iter_+1) ){
        init_alpha_value.clear();
    }
}

void JRCSAOPT::alpha_operation(int i)
{
    if( iter_count_ <= max_init_iter_ )alpha_operation_a(i);
    else alpha_operation_b(i);
}

void JRCSAOPT::alpha_operation_a(int i)
{
    arma::mat& alpha = *alpha_ptrlst_[i];
    alpha %= init_alpha_value[i];
    arma::vec alpha_rowsum = arma::sum(alpha,1);
    alpha.each_col() /= alpha_rowsum;
}

void JRCSAOPT::alpha_operation_b(int i)
{
    arma::uvec& input_label = *vll_ptrlst_[i];
    arma::mat& alpha = *alpha_ptrlst_[i];
    arma::uword lmax = input_label.max();
    arma::uword lmin = input_label.min();
    arma::fvec global_obj_prob_mean(obj_num_);
    std::vector<arma::uvec> oidxs(obj_num_);
    arma::rowvec global_alpha_mean = arma::mean(alpha);
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
        arma::mat current_alpha = alpha.rows(current_patch_index);
        arma::rowvec current_alpha_mean = arma::mean( current_alpha );
        arma::fvec current_obj_prob_mean( obj_num_ );
        #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            current_obj_prob_mean(o) = arma::accu( current_alpha_mean(oidxs[o]) );
            if( 0 != global_obj_prob_mean(o) )current_obj_prob_mean(o) /= global_obj_prob_mean(o);
        }
        arma::uword min_obj;
        assert( current_obj_prob_mean.is_finite() );
        current_obj_prob_mean.min(min_obj);
        assert(min_obj>=0&&min_obj<obj_num_);
        #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            if(o==min_obj)alpha(current_patch_index,oidxs[o]).fill(0.0);
            else alpha(current_patch_index,oidxs[o])*=current_obj_prob_mean(o);
        }
    }
    //normalise alpha
    arma::vec alpha_rowsum = arma::sum(alpha,1);
    alpha.each_col() /= alpha_rowsum;
}

void JRCSAOPT::rand_sphere(
        arma::fmat& ov
        )
{
    int k = ov.n_cols;
    arma::frowvec u = arma::randu<arma::frowvec>(k);
    arma::frowvec v = arma::randu<arma::frowvec>(k);
    arma::frowvec sign = arma::sign(arma::randu<arma::frowvec>(k) - 0.5);
    arma::frowvec cos_theta = 2*u-1;
    arma::frowvec phi = M_PI*(2*v-1);
    ov.row(0) = cos_theta%arma::cos(phi);
    ov.row(1) = cos_theta%arma::sin(phi);
    ov.row(2) = sign%arma::sqrt( ( 1.0 - arma::square(cos_theta) ) );
    #pragma omp parallel for
    for(arma::uword i=0;i<ov.n_cols;++i)
    {
        arma::fvec x = ov.col(i);
        ov.col(i) = x(randperm(3));
    }
}

}
