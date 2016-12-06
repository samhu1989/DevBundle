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
    connect(ui->UpdateCentroid,SIGNAL(clicked(bool)),this,SLOT(updateCentroid()));
    connect(ui->LoadFunc,SIGNAL(clicked(bool)),this,SLOT(loadFunc()));
    connect(ui->resEnergy,SIGNAL(clicked(bool)),this,SLOT(evaluateResEnergy()));
    connect(ui->SaveAll,SIGNAL(clicked(bool)),this,SLOT(saveAll()));

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
    if(re_func_.size()<inputs_.size())re_func_.resize(inputs_.size());
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
            if( s > e )
            {
                indices.push_back( sorted_i(i) );
                ++cnt;
            }else
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

void Spectrum::updateCentroid(void)
{
    arma::uword e = ui->MaxEigIndex->value();
    arma::uword index = 0;
    for(MeshList::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr m = *iter;
        arma::vec func = arma::min(arma::abs((*bases_[index]).cols(0,e)),1);
        ColorArray::colorfromValue((uint32_t*)m->custom_color_.vertex_colors(),func.size(),func);
        ++index;
    }
}

void Spectrum::evaluateResEnergy(void)
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
        arma::mat funcs;
        if(fname.endsWith(".mat"))
        {
            if(!MATIO::load_to_arma(funcs,fname.toStdString(),"X"))
            {
                std::cerr<<"Failed to load:"<<fname.toStdString()<<std::endl;
            }
        }else{
            funcs.load( fname.toStdString() );
        }
        arma::vec res(funcs.n_cols,arma::fill::zeros);
        for(int i=0;i<funcs.n_cols;++i)
        {
            res(i) = resEnergy(inputs_[index],index,funcs.col(i));
        }
        QFileInfo info(fname);
        QString res_name = info.path()+tr("/")+info.baseName()+tr("_res.mat");
        MATIO::save_to_matlab(res,res_name.toStdString(),"X");
        ++index;
    }
}

double Spectrum::resEnergy(MeshBundle<DefaultMesh>::Ptr m,const arma::uword index,const arma::vec& pixFunc)
{
    assert( m->mesh_.n_vertices() == pixFunc.size() );
    arma::uword s = ui->MinEigIndex->value();
    arma::uword e = ui->MaxEigIndex->value();
    arma::vec voxFunc;
    m->graph_.getVoxFunc(pixFunc,voxFunc);
    arma::vec coeff = arma::vectorise(voxFunc.t()*(*bases_[index]));
    arma::uvec sorted_i = arma::sort_index(arma::abs(coeff),"descend");
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
        if( s > e )//s > e then certianly and only use largest
        {
            indices.push_back( sorted_i(i) );
            ++cnt;
        }else
        if( sorted_i(i) > e || sorted_i(i) < s )
        {
            indices.push_back( sorted_i(i) );
            ++cnt;
        }
    }
    arma::uvec components(indices);
    arma::vec sub_coeff = coeff.rows(components);
    arma::mat sub_bases = (*bases_[index]).cols(components);
    voxFunc = sub_bases*sub_coeff;
    arma::vec tmp;
    m->graph_.getPixFunc(voxFunc,tmp);
    return arma::accu(arma::square(pixFunc - tmp));
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
    if(!re_func_[index])re_func_[index].reset(new arma::vec(m->mesh_.n_vertices(),arma::fill::randu));
    m->graph_.getPixFunc(voxFunc,*re_func_[index]);
    if(!res_[index])res_[index].reset(new arma::vec(m->mesh_.n_vertices(),arma::fill::randu));
    (*res_[index]) = (*func_[index]) - *re_func_[index];
    assert(re_func_[index]->is_finite());
    ColorArray::colorfromValue((uint32_t*)m->custom_color_.vertex_colors(),re_func_[index]->size(),*re_func_[index]);
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

void Spectrum::saveCoeff(QString fname)
{
    if(fname.isEmpty())
    {
        fname = QFileDialog::getSaveFileName(
                this,
                tr("Save Coefficients"),
                tr("../Dev_Data/"),
                tr(
                "Matlab Files (*.mat);;"
                "Armadillo Files (*.arma);;"
                "All Files (*)")
                );
    }
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

void Spectrum::saveLambda(QString fname)
{
    if(fname.isEmpty())
    {
        QString fname = QFileDialog::getSaveFileName(
                this,
                tr("Save Lambda"),
                tr("../Dev_Data/"),
                tr(
                "Matlab Files (*.mat);;"
                "Armadillo Files (*.arma);;"
                "All Files (*)")
                );
    }
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

void Spectrum::saveResidue(QString fname)
{
    if(fname.isEmpty())
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
    }
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

void Spectrum::saveReconstructed(QString fname)
{
    if(fname.isEmpty())
    {
        fname = QFileDialog::getSaveFileName(
                this,
                tr("Save Coefficients"),
                tr("../Dev_Data/"),
                tr(
                "Matlab Files (*.mat);;"
                "Armadillo Files (*.arma);;"
                "All Files (*)")
                );
    }
    if(fname.isEmpty())return;
    arma::mat func(re_func_[0]->size(),re_func_.size());
    arma::uword index = 0;
    for(std::vector<std::shared_ptr<arma::vec>>::iterator iter=re_func_.begin();iter!=re_func_.end();++iter)
    {
        func.col(index) = **iter;
        ++index;
    }
    if(fname.endsWith(".mat"))
    {
        MATIO::save_to_matlab<arma::mat>(func,fname.toStdString(),"X");
    }else{
        func.save( fname.toStdString() );
    }
}

void Spectrum::saveBases(QString fname)
{
    if(fname.isEmpty())
    {
        fname = QFileDialog::getSaveFileName(
                this,
                tr("Save Coefficients"),
                tr("../Dev_Data/"),
                tr(
                "Matlab Files (*.mat);;"
                "Armadillo Files (*.arma);;"
                "All Files (*)")
                );
    }
    if(fname.isEmpty())return;
    std::vector<std::string> name_lst_;
    name_lst_.reserve(inputs_.size());
    for(MeshList::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr m_ptr = *iter;
        name_lst_.push_back("base_"+m_ptr->name_);
    }
    if(fname.endsWith(".mat")){
        MATIO::save_to_matlab(bases_,fname.toStdString(),name_lst_);
    }else{
        MATIO::save_to_matlab(bases_,(fname+tr(".mat")).toStdString(),name_lst_);
    }
}

void Spectrum::saveAll(void)
{
    QString dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Save Coefficients"),
                tr("../Dev_Data/")
                );
    if(dirName.isEmpty())return;
    QString coeff_name,lambda_name,res_name,refunc_name,base_name;
    QDir dir;
    dir.setPath(dirName);
    coeff_name = dir.absoluteFilePath("./coeff.mat");
    lambda_name = dir.absoluteFilePath("./lambda.mat");
    res_name = dir.absoluteFilePath("./res.mat");
    refunc_name = dir.absoluteFilePath("./reconstructed.mat");
    base_name = dir.absoluteFilePath("./base.mat");
    saveCoeff(coeff_name);
    saveLambda(lambda_name);
    saveResidue(res_name);
    saveReconstructed(refunc_name);
    saveBases(base_name);
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
