#include "objectview.h"
#include "ui_objectview.h"
#include <armadillo>
ObjectView::ObjectView(
        std::vector<ObjModel::Ptr>& objects,
        QWidget *parent
        ):
    QFrame(parent),
    ui(new Ui::ObjectView),
    objects_(objects)
{
    ui->setupUi(this);
    widget_ = new MeshListViewerWidget();
    ui->gridLayout->addWidget(widget_);
    ui->color->setChecked(true);
    connect(ui->spatial_w,SIGNAL(toggled(bool)),this,SLOT(view_spatial_weight(bool)));
    connect(ui->color_w,SIGNAL(toggled(bool)),this,SLOT(view_color_weight(bool)));
    connect(ui->color,SIGNAL(toggled(bool)),this,SLOT(view_original_color(bool)));
    std::vector<ObjModel::Ptr>::iterator oiter;
    for(oiter=objects_.begin();oiter!=objects_.end();++oiter)
    {
        ObjModel& model = **oiter;
        widget_->list().push_back(model.GeoM_);
    }
}

void ObjectView::view_color_weight(bool view)
{
    if(view)
    {
        ui->color->setChecked(false);
        ui->color_w->setChecked(false);
        std::vector<ObjModel::Ptr>::iterator oiter;
        for(oiter=objects_.begin();oiter!=objects_.end();++oiter)
        {
            ObjModel& model = **oiter;
            arma::Col<uint32_t> var_color(
                        (uint32_t*)model.GeoM_->custom_color_.vertex_colors(),
                        model.GeoM_->mesh_.n_vertices(),
                        false,
                        true
                        );
            arma::fvec var(model.ColorP_);
//            arma::fvec geovar(model.GeoP_);
//            arma::uvec indices = arma::find(geovar==var);
//            std::cerr<<"equal num"<<indices.size()<<"/"<<var.size()<<std::endl;
            float max_var = arma::max(var);
            float min_var = arma::min(var);
            float h;
            ColorArray::RGB32 tmp;
            for(size_t idx=0;idx<var_color.size();++idx)
            {
                if(max_var!=min_var)h = ( var(idx) - min_var ) / ( max_var - min_var );
                else h = 0.0;
                if(!std::isfinite(h))std::logic_error("!std::isfinite(h)");
                if(h>1.0)std::logic_error("h>1.0");
                ColorArray::hsv2rgb(h*220.0+5.0,0.5,1.0,tmp);
                var_color(idx) = tmp.color;
            }
        }
        widget_->use_custom(true);
        widget_->updateGL();
    }
}

void ObjectView::view_spatial_weight(bool view)
{
    if(view)
    {
        ui->color->setChecked(false);
        ui->color_w->setChecked(false);
        std::vector<ObjModel::Ptr>::iterator oiter;
        for(oiter=objects_.begin();oiter!=objects_.end();++oiter)
        {
            ObjModel& model = **oiter;
            arma::Col<uint32_t> var_color(
                        (uint32_t*)model.GeoM_->custom_color_.vertex_colors(),
                        model.GeoM_->mesh_.n_vertices(),
                        false,
                        true
                        );
            arma::fvec var(model.GeoP_);
            float max_var = arma::max(var);
            float min_var = arma::min(var);
            float h;
            ColorArray::RGB32 tmp;
            for(size_t idx=0;idx<var_color.size();++idx)
            {
                if(max_var!=min_var)h = ( var(idx) - min_var ) / ( max_var - min_var );
                else h = 0.0;
                if(!std::isfinite(h))std::logic_error("!std::isfinite(h)");
                if(h>1.0)std::logic_error("h>1.0");
                ColorArray::hsv2rgb(h*220.0+5.0,0.5,1.0,tmp);
                var_color(idx) = tmp.color;
            }
        }
        widget_->use_custom(true);
        widget_->updateGL();
    }
}

void ObjectView::view_original_color(bool view)
{
    if(view)
    {
        widget_->use_custom(false);
        widget_->updateGL();
    }
}

ObjectView::~ObjectView()
{
    delete ui;
    widget_->deleteLater();
}
