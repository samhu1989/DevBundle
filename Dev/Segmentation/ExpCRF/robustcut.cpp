#include "robustcut.h"
#include <QMessageBox>
arma::umat RobustCut::base_segment_;
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

