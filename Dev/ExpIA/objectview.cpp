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
    connect(ui->norm_w,SIGNAL(toggled(bool)),this,SLOT(view_normal_weight(bool)));
    connect(ui->color,SIGNAL(toggled(bool)),this,SLOT(view_original_color(bool)));
    std::vector<ObjModel::Ptr>::iterator oiter;
    for(oiter=objects_.begin();oiter!=objects_.end();++oiter)
    {
        ObjModel& model = **oiter;
        widget_->list().push_back(model.GeoM_);
    }
    widget_->set_center_at_mesh(widget_->list().back()->mesh_);
    max_ = arma::fvec(objects_.size());
    min_ = arma::fvec(objects_.size());
    connect(&timer_,SIGNAL(timeout()),this,SLOT(update_labels()));
    timer_.setSingleShot(false);
    timer_.start(500);
    value_bar_ = QImage(ui->value_label->width(),ui->value_label->height(),QImage::Format_RGB888);
    for(size_t r=0;r<ui->value_label->height();++r)
    for(size_t c=0;c<ui->value_label->width();++c)
    {
        float h;
        h = float( ui->value_label->height() - r ) / float(ui->value_label->height());
        ColorArray::RGB32 tmp;
        ColorArray::hsv2rgb(h*220.0+5.0,0.5,1.0,tmp);
        value_bar_.setPixel(c,r,qRgb(tmp.rgba.r,tmp.rgba.g,tmp.rgba.b));
    }
    ui->value_label->setPixmap(QPixmap::fromImage(value_bar_));
}

void ObjectView::view_color_weight(bool view)
{
    if(view)
    {
        ui->color->setChecked(false);
        ui->color_w->setChecked(false);
        ui->norm_w->setChecked(false);
        std::vector<ObjModel::Ptr>::iterator oiter;
        size_t index = 0;
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
            max_(index) = max_var;
            min_(index) = min_var;
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
            ++index;
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
        ui->norm_w->setChecked(false);
        std::vector<ObjModel::Ptr>::iterator oiter;
        size_t index = 0;
        for(oiter=objects_.begin();oiter!=objects_.end();++oiter)
        {
            ObjModel& model = **oiter;
            arma::Col<uint32_t> var_color(
                        (uint32_t*)model.GeoM_->custom_color_.vertex_colors(),
                        model.GeoM_->mesh_.n_vertices(),
                        false,
                        true
                        );
            arma::fvec var(model.DistP_);
            float max_var = arma::max(var);
            float min_var = arma::min(var);
            max_(index) = max_var;
            min_(index) = min_var;
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
            ++index;
        }
        widget_->use_custom(true);
        widget_->updateGL();
    }
}

void ObjectView::view_normal_weight(bool view)
{
    if(view)
    {
        ui->color->setChecked(false);
        ui->color_w->setChecked(false);
        ui->spatial_w->setChecked(false);
        std::vector<ObjModel::Ptr>::iterator oiter;
        size_t index = 0;
        for(oiter=objects_.begin();oiter!=objects_.end();++oiter)
        {
            ObjModel& model = **oiter;
            arma::Col<uint32_t> var_color(
                        (uint32_t*)model.GeoM_->custom_color_.vertex_colors(),
                        model.GeoM_->mesh_.n_vertices(),
                        false,
                        true
                        );
            arma::fvec var(model.NormP_);
            float max_var = arma::max(var);
            float min_var = arma::min(var);
            max_(index) = max_var;
            min_(index) = min_var;
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
            ++index;
        }
        widget_->use_custom(true);
        widget_->updateGL();
    }
}

void ObjectView::view_original_color(bool view)
{
    if(view)
    {
        max_.fill(1.0);
        min_.fill(0.0);
        ui->spatial_w->setChecked(false);
        ui->color_w->setChecked(false);
        ui->norm_w->setChecked(false);
        widget_->use_custom(false);
        widget_->updateGL();
    }
}

void ObjectView::update_labels(void)
{
    uint32_t l = widget_->current_mesh_start();
    QString smin,smax;
    smin = smin.sprintf("%lf",min_(l));
    ui->min_label->setText(smin);
    smax = smax.sprintf("%lf",max_(l));
    ui->max_label->setText(smax);
}

ObjectView::~ObjectView()
{
    delete ui;
    widget_->deleteLater();
}
