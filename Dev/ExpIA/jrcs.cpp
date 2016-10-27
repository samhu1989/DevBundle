#include "jrcs.h"
#include "jrcsinitexternal.h"
#include "jrcsbase.h"
#include "jrcsaoni.h"
#include "jrcsaopt.h"
#include <QMessageBox>
#include "featurecore.h"
#include "iocore.h"
JRCSWork::DMatPtrLst JRCSWork::alpha_ptrlst_;
arma::fvec JRCSWork::obj_prob_;
JRCSWork::JRCSWork(
        MeshList& inputs,
        LabelList& labels,
        ModelList& objects,
        QObject *parent
        ) :
    inputs_(inputs),labels_(labels),objects_(objects),
    QObject(parent)
{
    ;
}

bool JRCSWork::configure(Config::Ptr config)
{
    if(labels_.empty())
    {
        QString msg = "Load a Label First\n";
        emit message(tr("JRCSWork Init:")+msg,0);
        return false;
    }
    if(inputs_.empty())
    {
        QString msg = "Load Input First\n";
        emit message(tr("JRCSWork Init:")+msg,0);
        return false;
    }
    if(inputs_[0]->graph_.empty())
    {
        QString msg = "Load Supervoxel First\n";
        emit message(tr("JRCSWork Init:")+msg,0);
        return false;
    }
    if(config->has("JRCS_obj_w"))
    {
        std::vector<float> obj_prob;
        config->getFloatVec("JRCS_obj_w",obj_prob);
        obj_prob_ = arma::fvec(obj_prob);
    }
    return true;
}

void JRCSWork::Init_SI_HKS()
{
    //extract HKS
    Feature::HKS<DefaultMesh>::MatPtrLst fLst;
    Feature::HKS<DefaultMesh> getFeature;
    fLst.resize(inputs_.size());
    Feature::HKS<DefaultMesh>::MatPtrLst::iterator fiter=fLst.begin();
    size_t index = 1;
    for(MeshBundle<DefaultMesh>::PtrList::const_iterator iter = inputs_.begin() ; iter != inputs_.end() ; ++iter )
    {
        QString path;
        path = path.sprintf("hks%02u.mat",index);
        (*fiter).reset(new arma::mat());
        emit message(path,0);
        getFeature.extract(*iter,**fiter);
        if(*fiter)MATIO::save_to_matlab(**fiter,(tr("./debug/HKS/")+path).toStdString(),"X");
        ++ fiter;
        ++ index;
        if( fiter == fLst.end() )break;
    }
    //extract the BOF of feature
    //    Feature::BOF bof;
    //    Feature::BOF::MatPtrLst histLst;
    //    bof.learn(fLst,labels_,histLst);
    //clustering on the BOF

    //generating obj_prob

    //generating alpha
    emit end();
}

void JRCSWork::Init_Bernolli()
{
    arma::uvec ref_label;
    Init_Bernolli_a(ref_label);
    Init_Bernolli_b(ref_label);
    emit end();
}

void JRCSWork::Init_Bernolli_a(arma::uvec& ref_label)
{
    std::cerr<<"Init_Bernolli_a"<<std::endl;
    arma::uvec sorted_i = arma::sort_index(obj_prob_);
    arma::fvec sort_obj_prob = obj_prob_(sorted_i);
    ref_label = arma::uvec(obj_prob_.size());
    obj_prob_ = obj_prob_(sorted_i);
    double min_dist_to_ref = std::numeric_limits<double>::max();
//    std::cerr<<"Find the existing segments that is closest to the specified object probability"<<std::endl;
    for(LabelList::iterator iter=labels_.begin();iter!=labels_.end();++iter)
    {
//        std::cerr<<"extract segmentation probability"<<std::endl;
        arma::uword min = arma::min(*iter);
        arma::uword max = arma::max(*iter);
        std::vector<float> seg_prob_vec;
        seg_prob_vec.reserve(max);
        std::vector<arma::uword> seg_label_vec;
        seg_label_vec.reserve(max);
        arma::uword label = min;
        while(label<=max)
        {
            arma::uvec indices = arma::find((*iter)==label);
            if(!indices.empty())
            {
                seg_label_vec.push_back(label);
                seg_prob_vec.push_back(float(indices.size())/float(iter->size()));
            }
            ++label;
        }
//        std:cerr<<"calculate the min_dist_to_ref and choose the ref"<<std::endl;
        arma::fvec seg_prob(seg_prob_vec);
        arma::uvec seg_label(seg_label_vec);
        arma::uvec sorted_j = arma::sort_index(seg_prob);
        arma::fvec sort_seg_prob = seg_prob(sorted_j);
        arma::uvec match_to_ref(ref_label.size());
        arma::fvec dists(sort_obj_prob);
        arma::sword j = sort_seg_prob.size() - 1;
//        std::cerr<<"a"<<std::endl;
        for( arma::sword i = sort_obj_prob.size() - 1 ; i >= 0  ; --i )
        {
//            std::cerr<<"j:"<<j<<std::endl;
//            std::cerr<<"i:"<<i<<std::endl;
            double new_dist = std::abs( sort_seg_prob(j) - sort_obj_prob(i)  );
            double old_dist = new_dist + std::numeric_limits<float>::epsilon();
            while( old_dist > new_dist )
            {
                dists(i) = new_dist;
//                std::cerr<<"j:"<<j<<std::endl;
//                std::cerr<<"seg_label.size():"<<seg_label.size()<<std::endl;
//                std::cerr<<"sorted_j("<<j<<"):"<<sorted_j(j)<<std::endl;
//                std::cerr<<"match_to_ref.size():"<<match_to_ref.size()<<std::endl;
                match_to_ref(i) = seg_label( sorted_j(j) );
                --j;
                if( j < 0 )break;
                old_dist = new_dist;
                new_dist = std::abs( sort_seg_prob(j) - sort_obj_prob(i) );
            }
        }
//        std::cerr<<"b"<<std::endl;
        if( arma::accu(dists) < min_dist_to_ref )
        {
            ref_label = match_to_ref;
        }
    }
}

void JRCSWork::Init_Bernolli_b(arma::uvec& ref_label)
{
    std::cerr<<"Init_Bernolli_b"<<std::endl;
    std::cerr<<"ref_label:"<<ref_label<<std::endl;
    alpha_ptrlst_.resize(inputs_.size());
    size_t index;
    MeshList::iterator iter;
    arma::uvec k_lst(inputs_.size());
    int idx = 0;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        k_lst(idx) = (*iter)->mesh_.n_vertices();
        ++idx;
    }
    arma::uword k_ = JRCS::JRCSBase::evaluate_k(k_lst);
    arma::uword r_k = k_;
    arma::uvec obj_size(obj_prob_.size());
    arma::fvec prob = obj_prob_;
    //determine the gaussian numbers for each object
    for(int oi = 0 ; oi < obj_prob_.size() ; ++ oi )
    {
        obj_size(oi) = arma::uword(float(k_)*float(prob(oi)));
        obj_size(oi) = std::max(arma::uword(9),obj_size(oi));
        obj_size(oi) = std::min(r_k,obj_size(oi));
        r_k -= obj_size(oi);
    }
    //calculate the column index by gaussian numbers
    std::vector<arma::uvec> cols;
    arma::uword sum = 0;
    for(int oi = 0 ; oi < obj_size.size() ; ++oi )
    {
        cols.emplace_back(obj_size[oi],arma::fill::zeros);
        cols.back() = arma::linspace<arma::uvec>( sum , sum + obj_size[oi] - 1 , obj_size[oi] );
        sum += obj_size[oi];
    }
    //fill alpha
    index = 0;
    for(DMatPtrLst::iterator iter = alpha_ptrlst_.begin() ; iter != alpha_ptrlst_.end() ; ++iter )
    {
        iter->reset(new arma::mat(inputs_[index]->mesh_.n_vertices(),k_,arma::fill::zeros));
        arma::mat& alpha = **iter;
        arma::uvec& label = labels_[index];
        alpha.fill(1.0/double(k_));
        #pragma omp parallel for
        for(int r=0;r<prob.n_rows;++r)
        {
            arma::uvec rows = arma::find( label == ref_label(r) );
            alpha(rows,cols[r]).fill(1.0);
        }
        ++index;
    }
}

bool JRCSWork::init_optimize(JRCSView* w)
{
    if( alpha_ptrlst_.empty() || obj_prob_.empty() )return false;
    std::shared_ptr<JRCS::JRCSInitBase> init_ptr(new JRCS::JRCSInitExternal(alpha_ptrlst_,obj_prob_));
    w->set_init_method(init_ptr);
    return true;
}

void JRCSWork::set_opt_aoni(JRCSView* w)
{
    std::shared_ptr<JRCS::JRCSBase> method(new JRCS::JRCSAONI());
    w->set_method(method);
}

void JRCSWork::set_opt_aopt(JRCSView* w)
{
    std::shared_ptr<JRCS::JRCSBase> method(new JRCS::JRCSAOPT());
    w->set_method(method);
}
