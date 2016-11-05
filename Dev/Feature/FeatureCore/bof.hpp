#include "bof.h"
namespace Feature {
BOF::BOF()
{
    codebook_size_ = 128;
}

void BOF::extract(const arma::mat& f,const arma::uvec& l,arma::mat& h)
{
    arma::uword label_max = arma::max(l);
    arma::uword label_min = arma::min(l);
    h = arma::mat(codebook_size_, label_max - label_min + 1 );
    //calculate idf
    arma::vec idf(gmm_.n_gaus(),arma::fill::zeros);
    arma::urowvec r = gmm_.assign(f,arma::eucl_dist);
    arma::mat counts(gmm_.n_gaus(),label_max,arma::fill::zeros);
    idf = arma::vec(gmm_.n_gaus(),arma::fill::zeros);
    for( arma::uword i=0 ; i < r.n_cols ; ++i )
    {
        if(l(i)==0)continue;
        assert( l(i) <= counts.n_cols );
        assert(r(i)<counts.n_rows);
        if( counts(r(i),l(i) - 1) == 0 ) counts( r(i) , l(i) - 1 ) = 1.0 ;
    }
    idf = arma::sum(counts,1);
    idf += 1.0;
    idf = label_max / idf;
    idf = arma::log(idf);
    //calculate tf-idf for each patch
    h.resize(gmm_.n_gaus(),label_max);
    arma::mat& tf = h;
    arma::rowvec word_num(label_max,arma::fill::ones);
    for( arma::uword i=0 ; i < r.n_cols ; ++i )
    {
        if(l(i)==0)continue;
        assert( l(i) <= tf.n_cols );
        assert( r(i) < tf.n_rows );
        tf( r(i) , l(i) - 1 ) += 1.0 ;
        word_num(l(i) - 1) += 1.0;
    }
    tf.each_row()/=word_num;
    tf.each_col()%=idf;
}

void BOF::learn(const MatPtrLst& f,const LabelLst& l,MatPtrLst& h)
{
    //learn code book by k-means
    arma::uword num = 0;
    arma::uword dim;
    arma::uvec start(f.size(),arma::fill::zeros);
    arma::uword count = start.size() - 1;
    for(MatPtrLst::const_iterator iter=f.cbegin() ; iter!=f.cend() ; ++iter )
    {
        dim = (*iter)->n_rows;
        start.tail(count) += (*iter)->n_cols;
        -- count;
        num += (*iter)->n_cols;
    }
    arma::mat data(dim,num);
    #pragma omp parallel for
    for(arma::uword i=0 ; i < f.size() ; ++i )
    {
        data.cols(start(i),start(i)+f[i]->n_cols - 1) = *f[i];
    }
    gmm_.learn(data,codebook_size_,arma::eucl_dist,arma::random_subset,50,0,1e-12,true);
    //calculate idf
    arma::uword index = 0;
    LabelLst::const_iterator liter = l.cbegin();
    idf_.resize(f.size());
    for(MatPtrLst::const_iterator iter = f.cbegin() ; iter != f.cend() ; ++iter )
    {
        arma::uword label_max = arma::max(*liter);
        arma::urowvec r = gmm_.assign(**iter,arma::eucl_dist);
        arma::mat counts(gmm_.n_gaus(),label_max,arma::fill::zeros);
        idf_[index] = arma::vec(gmm_.n_gaus(),arma::fill::zeros);
        for( arma::uword i=0 ; i < r.n_cols ; ++i )
        {
            if((*liter)(i)==0)continue;
            assert( (*liter)(i) <= counts.n_cols );
            assert( r(i) < counts.n_rows );
            if( counts( r(i) ,(*liter)(i) - 1) == 0 ) counts( r(i) , (*liter)(i) - 1 ) = 1.0 ;
        }
        idf_[index] = arma::sum(counts,1);
        idf_[index] += 1.0;
        idf_[index] = ( label_max + 1.0 ) / idf_[index];
        idf_[index] = arma::log( idf_[index] );
        ++index;
    }
    //calculate tf-idf for each patch
    h.resize(l.size());
    liter = l.cbegin();
    index = 0;
    for(MatPtrLst::const_iterator iter = f.cbegin() ; iter != f.cend() ; ++iter )
    {
        arma::uword label_max = arma::max(*liter);
        h[index].reset(new arma::mat(gmm_.n_gaus(),label_max,arma::fill::zeros));
        arma::urowvec r = gmm_.assign(**iter,arma::eucl_dist);
        arma::mat& tf = (*h[index]);
        arma::rowvec word_num(label_max,arma::fill::ones);
        for( arma::uword i=0 ; i < r.n_cols ; ++i )
        {
            if((*liter)(i)==0)continue;
            assert( (*liter)(i) <= tf.n_cols );
            assert( r(i) < tf.n_rows );
            tf( r(i) , (*liter)(i) - 1 ) += 1.0 ;
            word_num((*liter)(i) - 1) += 1.0;
        }
        tf.each_row()/=word_num;
        tf.each_col()%=idf_[index];
        ++index;
    }
}

}
