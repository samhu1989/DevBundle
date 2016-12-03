#include "spectrum.h"
#include "ui_spectrum.h"
#include <QFileDialog>
#include "iocore.h"
Spectrum::Spectrum(MeshList &inputs,Config::Ptr config,QWidget *parent) :
    inputs_(inputs),
    config_(config),
    QFrame(parent),
    ui(new Ui::Spectrum)
{
    ui->setupUi(this);
    connect(ui->UpdateSpec,SIGNAL(clicked(bool)),this,SLOT(updateSpectrum()));
    connect(ui->UpdateFunc,SIGNAL(clicked(bool)),this,SLOT(updateFunction()));
    connect(ui->LoadFunc,SIGNAL(clicked(bool)),this,SLOT(loadFunc()));
    connect(ui->SaveCoeff,SIGNAL(clicked(bool)),this,SLOT(saveCoeff()));
    connect(ui->SaveLambda,SIGNAL(clicked(bool)),this,SLOT(saveLambda()));
    connect(ui->SaveResidue,SIGNAL(clicked(bool)),this,SLOT(saveResidue()));

    ncuts_.configure(config_);
    ncuts_.setK( ui->EigNum->value() + 5 );
    ncuts_.setEpsilon(std::numeric_limits<double>::lowest());
    ncuts_.setType(Segmentation::NormalizedCuts<DefaultMesh>::M);

    updateBase();
    updateFunc();
}

void Spectrum::updateBase()
{
    if(bases_.size()<inputs_.size())bases_.resize(inputs_.size());
    if(coeff_.size()<inputs_.size())coeff_.resize(inputs_.size());
    if(func_.size()<inputs_.size())func_.resize(inputs_.size());
    if(lambda_.size()<inputs_.size())lambda_.resize(inputs_.size());
    if(res_.size()<inputs_.size())res_.resize(inputs_.size());
    arma::uword index = 0;
    ncuts_.setK( ui->EigNum->value() + 5 );
    for(MeshList::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr m_ptr = *iter;
        if( !bases_[index] || bases_[index]->n_cols != ui->EigNum->value()+1 )
        {
            ncuts_.computeW_Graph(m_ptr);
            ncuts_.decompose();
            bases_[index].reset(new arma::mat());
            (*bases_[index]) = ncuts_.getY().cols(0,ui->EigNum->value());
            lambda_[index].reset(new arma::vec());
            (*lambda_[index]) = ncuts_.getLambda().rows(0,ui->EigNum->value());
            assert((*bases_[index]).is_finite());
            assert((*lambda_[index]).is_finite());
        }
        ++index;
    }
}

void Spectrum::updateFunc(void)
{
    arma::uword s = ui->MinEigIndex->value();
    arma::uword e = ui->MaxEigIndex->value();
    if( s > e )e = s;
    arma::uword index = 0;
    for(MeshList::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr m_ptr = *iter;
        if(!func_[index])func_[index].reset(new arma::vec(m_ptr->mesh_.n_vertices(),arma::fill::randu));
        arma::vec pixFunc = *func_[index];
        arma::vec voxFunc;
        m_ptr->graph_.getVoxFunc(pixFunc,voxFunc);
        if(!coeff_[index])coeff_[index].reset(new arma::vec(ui->EigNum->value() + 1 ));
        (*coeff_[index]) = arma::vectorise(voxFunc.t()*(*bases_[index]));
        arma::vec coeff = arma::abs( *coeff_[index] );
        arma::uvec sorted_i = arma::sort_index(coeff,"descend");
        std::vector<arma::uword> indices;
        for( arma::uword i=0 ; i < sorted_i.size() ; ++i )
        {
            if( sorted_i(i) <= e && sorted_i(i) >= s )
            {
                indices.push_back( sorted_i(i) );
            }
        }
        arma::uword cnt = 0;
        for( arma::uword i=0 ; i < sorted_i.size() && cnt < ui->Additional->value(); ++i )
        {
            if( sorted_i(i) > e || sorted_i(i) < s )
            {
                indices.push_back( sorted_i(i) );
                ++cnt;
            }
        }
        composeFunc(m_ptr,index,arma::uvec(indices));
        ++ index;
    }
}

void Spectrum::composeFunc(
        const MeshBundle<DefaultMesh>::Ptr m,
        const arma::uword index,
        const arma::uvec& components
        )
{
    arma::vec sub_coeff = (*coeff_[index]).rows(components);
    arma::mat sub_bases = (*bases_[index]).cols(components);
    arma::vec voxFunc = sub_bases*sub_coeff;
    arma::vec pixFunc;
    m->graph_.getPixFunc(voxFunc,pixFunc);
    if(!res_[index])res_[index].reset(new arma::vec(m->mesh_.n_vertices(),arma::fill::randu));
    (*res_[index]) = (*func_[index]) - pixFunc;
    assert(pixFunc.is_finite());
    ColorArray::colorfromValue((uint32_t*)m->custom_color_.vertex_colors(),pixFunc.size(),pixFunc);
}

void Spectrum::loadFunc(void)
{
    QStringList fileNames = QFileDialog::getOpenFileNames(this,
        tr("Open function files"),
        tr("../Dev_Data/"),
        tr(
        "Matlab Files (*.mat);;"
        "Armadillo Files (*.arma);;"
        "All Files (*)")
    );
    if(fileNames.empty())return;
    arma::uword index = 0;
    foreach(QString fname,fileNames)
    {
        if( index >= func_.size() )break;
        if(fname.endsWith(".mat"))
        {
            if(!MATIO::load_to_arma<arma::vec>(*func_[index],fname.toStdString(),"X"))
            {
                std::cerr<<"Failed to load:"<<fname.toStdString()<<std::endl;
            }
        }else{
            func_[index]->load( fname.toStdString() );
        }
        ++index;
    }
}

void Spectrum::saveCoeff(void)
{
    QString fname = QFileDialog::getSaveFileName(
                this,
                tr("Save Coefficients"),
                tr("../Dev_Data/"),
                tr(
                "Matlab Files (*.mat);;"
                "Armadillo Files (*.arma);;"
                "All Files (*)")
                );
    if(fname.isEmpty())return;
    arma::mat coeffs(coeff_[0]->size(),coeff_.size());
    arma::uword index = 0;
    for(std::vector<std::shared_ptr<arma::vec>>::iterator iter=coeff_.begin();iter!=coeff_.end();++iter)
    {
        coeffs.col(index) = **iter;
        ++index;
    }
    if(fname.endsWith(".mat"))
    {
        MATIO::save_to_matlab<arma::mat>(coeffs,fname.toStdString(),"X");
    }else{
        coeffs.save( fname.toStdString() );
    }
}

void Spectrum::saveLambda(void)
{
    QString fname = QFileDialog::getSaveFileName(
                this,
                tr("Save Coefficients"),
                tr("../Dev_Data/"),
                tr(
                "Matlab Files (*.mat);;"
                "Armadillo Files (*.arma);;"
                "All Files (*)")
                );
    if(fname.isEmpty())return;
    arma::mat lambda(lambda_[0]->size(),lambda_.size());
    arma::uword index = 0;
    for(std::vector<std::shared_ptr<arma::vec>>::iterator iter=lambda_.begin();iter!=lambda_.end();++iter)
    {
        lambda.col(index) = **iter;
        ++index;
    }
    if(fname.endsWith(".mat"))
    {
        MATIO::save_to_matlab<arma::mat>(lambda,fname.toStdString(),"X");
    }else{
        lambda.save( fname.toStdString() );
    }
}

void Spectrum::saveResidue(void)
{
    QString fname = QFileDialog::getSaveFileName(
                this,
                tr("Save Coefficients"),
                tr("../Dev_Data/"),
                tr(
                "Matlab Files (*.mat);;"
                "Armadillo Files (*.arma);;"
                "All Files (*)")
                );
    if(fname.isEmpty())return;
    arma::mat res(res_[0]->size(),res_.size());
    arma::uword index = 0;
    for(std::vector<std::shared_ptr<arma::vec>>::iterator iter=res_.begin();iter!=res_.end();++iter)
    {
        res.col(index) = **iter;
        ++index;
    }
    if(fname.endsWith(".mat"))
    {
        MATIO::save_to_matlab<arma::mat>(res,fname.toStdString(),"X");
    }else{
        res.save( fname.toStdString() );
    }
}

void Spectrum::updateFunction(void)
{
    arma::uword index = 0;
    for(MeshList::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr m_ptr = *iter;
        if(!func_[index])func_[index].reset(new arma::vec(m_ptr->mesh_.n_vertices(),arma::fill::randu));
        arma::vec pixFunc = *func_[index];
        ColorArray::colorfromValue((uint32_t*)m_ptr->custom_color_.vertex_colors(),pixFunc.size(),pixFunc);
        ++ index;
    }
}

void Spectrum::updateSpectrum(void)
{
    updateBase();
    updateFunc();
}

Spectrum::~Spectrum()
{
    delete ui;
}
