#include "robustcut.h"
#include <QMessageBox>
#include <QPainter>
arma::umat RobustCut::base_segment_;
RobustCut::ColorLabelMap RobustCut::map_;
void RobustCut::base_to_image(const arma::uvec& label_,QImage& img_)
{
    for(int y=0;y<img_.height();++y)
        for(int x=0;x<img_.width();++x)
        {
            int h,s,l,ch_l;
            QColor cv(img_.pixel(x,y));
            if( map_.find( label_(x+y*img_.width()) )==map_.end() )
            {
                std::srand(label_(x+y*img_.width())+1);
                int index = std::rand()%( ColorArray::DefaultColorNum_ - 1 );
                index += 1;
                map_.insert(label_(x+y*img_.width()),arma::uword(ColorArray::DefaultColor[index].color));
            }
            arma::uword ci = map_.value( label_(x+y*img_.width()) );
            QColor ch = QColor::fromRgb(ci);
            cv.getHsv(&h,&s,&l);
            ch.getHsv(&h,&s,&ch_l);
            img_.setPixel(x,y,ch.rgb());
        }
}
void RobustCut::save_base_to_image(uint32_t w,uint32_t h,QImage& img)
{
    uint32_t n = base_segment_.n_cols;
    uint32_t nc = std::sqrt(n);
    uint32_t nr = n / nc;
    if(n>(nr*nc))nr+=1;
    img = QImage(w*nc+nc-1,h*nr+nr-1,QImage::Format_RGB888);
    img.fill(Qt::black);
    QImage tmp(w,h,QImage::Format_RGB888);
    QPainter paint;
    paint.begin(&img);
    for(arma::uword i=0;i<n;++i)
    {
        int ny = i / nc;
        int nx = i - ny*nc;
        int x = nx*(w+1);
        int y = ny*(h+1);
        base_to_image(base_segment_.col(i),tmp);
        paint.drawImage(QPoint(x,y),tmp);
    }
    paint.end();
}
RobustCut::RobustCut(
        QImage &inputs,
        arma::uvec &labels,
        QObject *parent
        ):inputs_(inputs),labels_(labels),QObject(parent)
{
    setObjectName("RobustCut");
    base_segment_i_ = 0;
    if(!base_segment_.empty())
    {
        base_segment_N_ = base_segment_.n_cols;
    }else base_segment_N_ = 0;
}

bool RobustCut::configure(Config::Ptr config)
{
    if(inputs_.isNull())return false;
    return cuts_.configure(config);
}

void RobustCut::base_segments()
{
    base_segment_.clear();
    base_segment_i_ = 0;
    QString msg;
    std::cerr<<"============"<<std::endl;
    cuts_.generate_base_segment(inputs_);
    std::cerr<<"============"<<std::endl;
    base_segment_ = cuts_.base_segments();
    base_segment_N_=base_segment_.n_cols;
    emit message(msg.sprintf("Done Base Segment"),0);
    show_base_segment();
}

void RobustCut::keyPressEvent(QKeyEvent* event)
{
    if(event->type()!=QEvent::KeyPress)return;
    if(event->key()==Qt::Key_Escape)exit_base_segment();
    if(event->key()==Qt::Key_Left)last_base_segment();
    if(event->key()==Qt::Key_Right)next_base_segment();
}

void RobustCut::next_base_segment()
{
    base_segment_i_ ++;
    if(base_segment_i_==base_segment_N_)base_segment_i_=0;
    show_base_segment();
}

void RobustCut::last_base_segment()
{
    if(base_segment_i_==0)base_segment_i_ = base_segment_N_ - 1;
    else base_segment_i_ --;
    show_base_segment();
}

void RobustCut::exit_base_segment()
{
    emit end();
}

void RobustCut::show_base_segment()
{
    assert(base_segment_i_<base_segment_.n_cols);
    labels_ = base_segment_.col(base_segment_i_);
    labels_ += 1;
    QString msg;
    emit message(msg.sprintf("Current:%u/%u,Next:->,Last:<-,Esc:Exit",base_segment_i_+1,base_segment_N_),0);
}

void RobustCut::consensus_segment()
{
    cuts_.base_segments() = base_segment_;
    cuts_.solve_consensus_segment(inputs_,labels_);
    emit end();
}

