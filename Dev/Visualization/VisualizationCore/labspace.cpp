#include "labspace.h"
#include "ui_labspace.h"
#include <QPainter>
#include "common.h"
#include <QSizePolicy>
LabSpace::LabSpace(QWidget *parent) :
    QFrame(parent),
    ui(new Ui::LabSpace)
{
    ui->setupUi(this);
    ab_ = new LabLabel(LabLabel::ab);
    L_ = new LabLabel(LabLabel::L);
    ui->horizontalLayout->addWidget(ab_);
    ui->horizontalLayout->addWidget(L_);
}

LabSpace::~LabSpace()
{
    std::cerr<<"s"<<std::endl;
    ui->horizontalLayout->removeWidget(ab_);
    ui->horizontalLayout->removeWidget(L_);
    std::cerr<<"s"<<std::endl;
    ab_->deleteLater();
    L_->deleteLater();
    std::cerr<<"s"<<std::endl;
    delete ui;
    std::cerr<<"s"<<std::endl;
}
