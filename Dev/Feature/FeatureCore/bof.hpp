#include "bof.h"
namespace Feature {
BOF::BOF()
{
    codebook_size_ = 128;
}

bool BOF::configure(Config::Ptr)
{
    ;
}

void BOF::extract(const arma::mat& f,const arma::uvec& l,arma::mat& h)
{
    arma::uword label_max = arma::max(l);
    arma::uword label_min = arma::min(l);
    h = arma::mat(codebook_size_, label_max - label_min + 1 );
    #pragma omp parallel for
    for(arma::uword label=label_min ; label <= label_max ; ++label )
    {
        arma::uvec indices = arma::find(l==label);
        if(!indices.empty())
        {
            h.col(label) = extract(f.cols(indices)) ;
        }else{
            h.col(label).fill(0.0);
        }
    }
}

arma::vec BOF::extract(const arma::mat& f)
{
    arma::vec result(gmm_.n_gaus(),arma::fill::zeros);
    arma::urowvec c = gmm_.assign(f,arma::eucl_dist);
    for(arma::urowvec::iterator iter = c.begin() ; iter != c.end() ; ++ iter)
    {
        result(*iter) += 1.0;
    }
    result /= c.size();
    result %= idf_;
    return result;
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
    idf_ = arma::vec(gmm_.n_gaus(),arma::fill::zeros);
    double doc_num = 0;
    LabelLst::const_iterator liter = l.cbegin();
    for(MatPtrLst::const_iterator iter = f.cbegin() ; iter != f.cend() ; ++iter )
    {
        arma::uword label_max = arma::max(*liter);
        arma::uword label_min = arma::min(*liter);
        arma::urowvec r = gmm_.assign(**iter,arma::eucl_dist);
        arma::mat counts(label_max - label_min + 1,gmm_.n_gaus());
        for( arma::uword i=0 ; i < r.n_cols ; ++i )
        {
            if( counts((*liter)(i),r(i)) == 0 ) counts((*liter)(i),r(i)) = 1.0 ;
        }
        doc_num += label_max - label_min + 1;
        idf_ += arma::sum(counts,1);
    }
    idf_ = doc_num / idf_;
    idf_ = arma::log(idf_);
}

}
