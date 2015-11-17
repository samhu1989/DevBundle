#ifndef LABSPACE_H
#define LABSPACE_H
#include <QFrame>
#include <QPaintEvent>
#include <armadillo>
#include "visualizationcore_global.h"
#include <QLabel>
#include <QPainter>
namespace Ui {
class LabSpace;
}
class LabLabel:public QLabel
{
    Q_OBJECT
public:
    typedef enum{
        L,ab
    }Mode;
    explicit LabLabel(Mode m,QWidget *parent = 0):
        QLabel(parent),
        m_(m)
    {
        switch(m_)
        {
        case L:
            setMinimumSize(50,200);
            setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Fixed);
            L_ = QImage(50,100,QImage::Format_RGB888);
            L_.fill(Qt::magenta);
            std::cerr<<"L size:"<<L_.byteCount()<<std::endl;
            PL_ = arma::Mat<uint8_t>((uint8_t*)L_.bits(),3,L_.byteCount()/3,false,true);
            break;
        case ab:
            setMinimumSize(200,200);
            setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Fixed);
            ab_ = QImage(256,256,QImage::Format_RGB888);
            std::cerr<<"ab size:"<<ab_.byteCount()<<std::endl;
            L50ab_ = arma::Mat<uint8_t>((uint8_t*)ab_.bits(),3,ab_.byteCount()/3,false,true);
            break;
        }
    }

    ~LabLabel()
    {
        std::cerr<<"-"<<std::endl;
    }

protected:
    void paintEvent(QPaintEvent *e)
    {
        QPainter p;
        p.begin(this);
        switch(m_)
        {
        case L:
            p.drawImage(0,0,L_.scaledToHeight(this->height()));
            break;
        case ab:
            p.drawImage(0,0,ab_.scaled(this->size()));
            break;
        }
        p.end();
        QLabel::paintEvent(e);
    }
private:
    QImage ab_;
    QImage L_;
    arma::Mat<uint8_t> L50ab_;
    arma::Mat<uint8_t> PL_;
    Mode m_;
};

class VISUALIZATIONCORESHARED_EXPORT LabSpace : public QFrame
{
    Q_OBJECT
public:
    explicit LabSpace(QWidget *parent = 0);
    ~LabSpace();
private:
    LabLabel* ab_;
    LabLabel* L_;
    Ui::LabSpace *ui;
};

#endif // LABSPACE_H
