#include "segview.h"
#include "ui_segview.h"
#include <QColor>
SegView::SegView(QImage& img,
        arma::uvec& lbl, ColorLabelMap &map,
        QWidget *parent
        ):
    img_(img),
    label_(lbl),
    map_(map),
    view_label_(false),
    QWidget(parent),
    ui(new Ui::SegView)
{
    ui->setupUi(this);
    setMinimumSize(320,240);
    setWindowTitle(img_.text(tr("Path")));
    ui->label->setAlignment(Qt::AlignHCenter|Qt::AlignCenter);
}

void SegView::paintEvent(QPaintEvent* e)
{
    if(view_label_&&label_.size()==(img_.height()*img_.width()))
    {
        QImage img(img_.width(),img_.height(),QImage::Format_RGB888);
        for(int y=0;y<img.height();++y)
            for(int x=0;x<img.width();++x)
            {
                int h,s,l,ch_l;
                QColor cv = img_.pixelColor(x,y);
                arma::uword ci = map_.key( label_(x+y*img_.width()) );
                QColor ch = QColor::fromRgb(ci);
                cv.getHsl(&h,&s,&l);
                ch.getHsl(&h,&s,&ch_l);
                img.setPixelColor(x,y,QColor::fromHsl(h,s,(ch_l+l)/2));
            }
        ui->label->setPixmap(QPixmap::fromImage(img));
    }else{
        ui->label->setPixmap(QPixmap::fromImage(img_));
    }
    QWidget::paintEvent(e);
}

SegView::~SegView()
{
    delete ui;
}
