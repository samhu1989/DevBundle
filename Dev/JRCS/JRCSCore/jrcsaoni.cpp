#include "jrcsaoni.h"
namespace JRCS {
JRCSAONI::JRCSAONI():JRCSBase()
{

}
void JRCSAONI::computeOnce()
{
    //reset sum
    xv_sum_.fill(0.0);
    xn_sum_.fill(0.0);
    xc_sum_.fill(0.0);
    var_sum.fill(0.0);
    alpha_sum.fill(0.0);
    alpha_sumij.fill(0.0);

    //reset transformed latent center
    xtc_ = *xc_ptr_;


    std::vector<arma::uvec> oidx(obj_num_);

    #pragma omp parallel for
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        oidx[o] = arma::find(obj_label_==(o+1));
    }

    for( int i = 0 ; i < vvs_ptrlst_.size() ; ++ i )
    {
        float rate0 = float(i) / float(i+1);
        float rate1 = 1.0 - rate0;

    //reset transformed latent center
        xtv_ = *xv_ptr_;
        xtn_ = *xn_ptr_;

    //input vertices normal color
        arma::fmat& vv_ = *vvs_ptrlst_[i];
        arma::fmat& vn_ = *vns_ptrlst_[i];
        arma::Mat<uint8_t>& vc_ = *vcs_ptrlst_[i];
    //update alpha
        arma::fmat& alpha = *alpha_ptrlst_[i];
        Ts& rt = rt_lst_[i];

    //transform the object (transform the latent center one object by one object)
        #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            arma::fmat R(rt[o].R,3,3,false,true);
            arma::fvec t(rt[o].t,3,false,true);
            arma::fmat& objv = *objv_ptrlst_[o];
            arma::fmat& objn = *objn_ptrlst_[o];
            objv = R*objv;
            objv.each_col() += t;
            objn = R*objn;
        }

//        arma::fmat tmpxc = arma::conv_to<arma::fmat>::from(xtc_);
//        arma::fmat tmpvc = arma::conv_to<arma::fmat>::from(vc_);

   //calculate alpha
        if( (iter_count_>0) || (!init_alpha_) )
        {
            #pragma omp parallel for
            for(int r = 0 ; r < alpha.n_rows ; ++r )
            {
                arma::fmat tmpv = xtv_.each_col() - vv_.col(r);
                alpha.row(r)  = arma::sum(arma::square(tmpv));
            }

            alpha.each_row() %= (-0.5*x_invvar_);
            alpha = arma::exp(alpha);
            alpha.each_row() %= arma::pow(x_invvar_,1.5);
            alpha.each_row() %= x_p_;

            #pragma omp parallel for
            for(int o = 0 ; o < obj_num_ ; ++o )
            {
                alpha.cols(oidx[o]) *= obj_prob_(o);
            }
        }else{
            if(verbose_>0)std::cerr<<"using init alpha"<<std::endl;
        }

    //normalise alpha
        arma::fvec alpha_rowsum = ( 1.0 + beta_ ) * arma::sum(alpha,1);
        assert(alpha_rowsum.is_finite());
        alpha.each_col() /= alpha_rowsum;
        alpha_rowsum = arma::sum(alpha,1);
        assert(alpha_rowsum.is_finite());
    //constraint alpha
        if(iter_count_>0)alpha_operation(i);

    //update RT
    //#1 calculate weighted point cloud
        arma::fmat& wv = *wvs_ptrlst_[i];
        arma::fmat& wn = *wns_ptrlst_[i];
        arma::fmat wc = arma::conv_to<arma::fmat>::from(*wcs_ptrlst_[i]);
        arma::frowvec alpha_colsum = arma::sum( alpha );
        arma::frowvec alpha_median = arma::median( alpha );
        arma::fmat trunc_alpha = alpha;
        assert(trunc_alpha.is_finite());

        #pragma omp parallel for
        for(int c=0;c<alpha.n_cols;++c)
        {
            arma::fvec col = trunc_alpha.col(c);
            col( col < alpha_median(c) ).fill(0.0);
            assert(col.is_finite());
            trunc_alpha.col(c) = col;
        }

        arma::frowvec trunc_alpha_colsum = arma::sum(trunc_alpha);
        assert(trunc_alpha_colsum.is_finite());

        wv = vv_*trunc_alpha;
        wn = vn_*trunc_alpha;
        wc = arma::conv_to<arma::fmat>::from(vc_)*trunc_alpha;

        #pragma omp parallel for
        for(int c=0;c<alpha.n_cols;++c)
        {
            if( 0 < trunc_alpha_colsum(c) )
            {
                wv.col(c) /= trunc_alpha_colsum(c);
                wn.col(c) /= trunc_alpha_colsum(c);
                wc.col(c) /= trunc_alpha_colsum(c);
            }
        }

        wn = arma::normalise( wn );
        *wcs_ptrlst_[i] = arma::conv_to<arma::Mat<uint8_t>>::from(wc);

        #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            arma::fmat A;
            arma::fmat U,V;
            arma::fvec s;
            arma::fmat R(rt[o].R,3,3,false,true);
            arma::fvec t(rt[o].t,3,false,true);
            arma::fmat dR;
            arma::fvec dt;
            arma::fmat objv = *objv_ptrlst_[o];
            arma::fmat v;
            v = wv.cols(oidx[o]);
            assert(v.is_finite());
            arma::fmat cv = v.each_col() - arma::mean(v,1);
            objv.each_col() -= arma::mean(objv,1);
            A = cv*objv.t();
            assert(A.is_finite());
            switch(rttype_)
            {
            case Gamma:
            {
                arma::fmat B = A.submat(0,0,1,1);
                dR = arma::fmat(3,3,arma::fill::eye);
                if(arma::svd(U,s,V,B,"std"))
                {
                    arma::fmat C(2,2,arma::fill::eye);
                    C(1,1) = arma::det( U * V.t() )>=0 ? 1.0 : -1.0;
                    arma::fmat dR2D = U*C*(V.t());
                    dR.submat(0,0,1,1) = dR2D;
                    dt = arma::mean( v - dR*(*objv_ptrlst_[o]),1);
                }
            }
                break;
            default:
            {
                if(arma::svd(U,s,V,A,"std"))
                {
                    arma::fmat C(3,3,arma::fill::eye);
                    C(2,2) = arma::det( U * V.t() )>=0 ? 1.0 : -1.0;
                    dR = U*C*(V.t());
                    dt = arma::mean( v - dR*(*objv_ptrlst_[o]),1);
                }
            }
            }

            //updating objv
            //not needed since the transformed latent center is reset at the begining for each frame
            /*
            *objv_ptrlst_[o] = dR*(*objv_ptrlst_[o]);
            (*objv_ptrlst_[o]).each_col() += dt;
            *objn_ptrlst_[o] = dR*(*objn_ptrlst_[o]);
            */

            assert(dR.is_finite());
            assert(dt.is_finite());

            //updating R T
            R = dR*R;
            t = dR*t + dt;

            //accumulate for updating X
            arma::fmat tv = v.each_col() - t;
            xv_sum_.cols(oidx[o]) =  rate0*xv_sum_.cols(oidx[o])+rate1*(R.i()*tv);
            xn_sum_.cols(oidx[o]) = rate0*xn_sum_.cols(oidx[o])+rate1*R.i()*wn.cols(oidx[o]);
            xc_sum_.cols(oidx[o]) = rate0*xc_sum_.cols(oidx[o])+rate1*wc.cols(oidx[o]);
            //            }
        }
        //update var
        alpha_sum += alpha_colsum;
        arma::fmat alpha_2(alpha.n_rows,alpha.n_cols);
        #pragma omp parallel for
        for(int r=0;r<alpha_2.n_rows;++r)
        {
            alpha_2.row(r) = arma::sum(arma::square(xtv_.each_col() - vv_.col(r)));
        }
        arma::frowvec tmpvar = arma::sum(alpha_2%alpha);
        var_sum += tmpvar;
        alpha_sumij += alpha_colsum;
        QCoreApplication::processEvents();
    }
    //Updating X
    assert(xv_sum_.is_finite());
    *xv_ptr_ = xv_sum_;
//    assert((*xv_ptr_).has_inf()||(*xv_ptr_).has_nan());
    //fix the x center position
    #pragma omp parallel for
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        arma::fmat newxv = xv_ptr_->cols(oidx[o]);
        arma::fvec t =  obj_pos_.col(o) - arma::mean(newxv,1);
        xv_ptr_->cols(oidx[o]) = newxv.each_col() + t;
    }

    *xn_ptr_ = xn_sum_;
    *xn_ptr_ = arma::normalise(*xn_ptr_);
    *xc_ptr_ = arma::conv_to<arma::Mat<uint8_t>>::from( xc_sum_ );
    //restore reciprocal fractions of variation
    x_invvar_ = ( (3.0*alpha_sum ) / ( var_sum + beta_ ) );

    //red means low in invvar and high in var
    //blue means opposite
    ColorArray::colorfromValue((ColorArray::RGB888*)xc_ptr_->memptr(),xc_ptr_->n_cols,x_invvar_.t());

    float mu = arma::accu(alpha_sumij);
    x_p_ = alpha_sumij;

    if( mu != 0)x_p_ /= mu;
    #pragma omp parallel for
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        obj_prob_(o) = arma::accu( x_p_(oidx[o]) );
    }
    obj_prob_ += beta_; //add a small number to prevent underflow
    obj_prob_ /= arma::accu(obj_prob_);
    if(verbose_>0)std::cerr<<"obj_prob:"<<obj_prob_<<std::endl;
}
}
