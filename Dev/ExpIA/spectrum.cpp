#include "spectrum.h"
#include "ui_spectrum.h"

Spectrum::Spectrum(MeshList &inputs,Config::Ptr config,QWidget *parent) :
    inputs_(inputs),
    config_(config),
    QFrame(parent),
    ui(new Ui::Spectrum)
{
    ui->setupUi(this);
    connect(ui->UpdateFunc,SIGNAL(clicked(bool)),this,SLOT(update()));

    ncuts_.configure(config_);
    ncuts_.setK(ui->EigNum->value());
    ncuts_.setType(Segmentation::NormalizedCuts<DefaultMesh>::M);

    updateBase();
    updateFunc();

}

void Spectrum::updateBase()
{
    if(bases_.size()<inputs_.size())bases_.resize(inputs_.size());
    if(coeff_.size()<inputs_.size())coeff_.resize(inputs_.size());
    if(func_.size()<inputs_.size())func_.resize(inputs_.size());
    arma::uword index = 0;
    for(MeshList::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr m_ptr = *iter;
        if( !bases_[index] || bases_[index]->n_cols < ( ui->EigNum->value() + 1 ) )
        {
            ncuts_.computeW_Graph(m_ptr);
            ncuts_.decompose();
            bases_[index].reset(new arma::mat());
            (*bases_[index]) = ncuts_.getY();
        }
        ++index;
    }
}

void Spectrum::updateFunc()
{
    arma::uword index = 0;
    arma::uword s = ui->MinEigIndex->value();
    arma::uword e = ui->MaxEigIndex->value();
    if( s > e )e = s;
    for(MeshList::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr m_ptr = *iter;
        if(!func_[index])func_[index].reset(new arma::vec(m_ptr->mesh_.n_vertices(),arma::fill::randu));
        arma::vec pixFunc = *func_[index];
        arma::vec voxFunc;
        toVoxFunc(m_ptr,pixFunc,voxFunc);
        if(!coeff_[index])coeff_[index].reset(new arma::vec(voxFunc.n_rows));
        (*coeff_[index]) = arma::vectorise(voxFunc.t()*(*bases_[index]));
        arma::vec sub_coeff = (*coeff_[index]).rows(s,e);
        arma::mat sub_bases = (*bases_[index]).cols(s,e);
        voxFunc = sub_bases*sub_coeff;
        toPixFunc(m_ptr,voxFunc,pixFunc);
        ColorArray::colorfromValue((uint32_t*)m_ptr->custom_color_.vertex_colors(),pixFunc.size(),pixFunc);
        ++ index;
    }
}

void Spectrum::loadFunc(void)
{
    ;
}

void Spectrum::loadCoeff(void)
{
    ;
}

void Spectrum::update(void)
{
    updateBase();
    updateFunc();
}

void Spectrum::toPixFunc(MeshBundle<DefaultMesh>::Ptr m,const arma::vec& voxFunc,arma::vec& pixFunc)
{
    pixFunc = arma::vec(m->mesh_.n_vertices(),arma::fill::zeros);
    arma::uvec indices = m->graph_.voxel_label;
    indices -= 1;
    pixFunc = voxFunc(indices);
}

void Spectrum::toVoxFunc(MeshBundle<DefaultMesh>::Ptr m,const arma::vec& pixFunc,arma::vec& voxFunc)
{
    voxFunc = arma::vec(m->graph_.size(),arma::fill::zeros);
    arma::vec voxCunt = arma::vec(m->graph_.size(),arma::fill::zeros);
    const arma::uvec& indices = m->graph_.voxel_label;
    arma::vec::const_iterator piter = pixFunc.cbegin();
    for(arma::uvec::const_iterator iiter=indices.cbegin();iiter!=indices.cend();++iiter)
    {
        voxFunc( (*iiter - 1) ) = pixFunc(*piter);
        voxCunt( (*iiter - 1) ) += 1.0;
        ++ piter;
        if( piter == pixFunc.cend() )break;
    }
    voxFunc /= voxCunt;
}

Spectrum::~Spectrum()
{
    delete ui;
}
