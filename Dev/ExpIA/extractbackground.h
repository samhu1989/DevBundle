#ifndef EXTRACTBACKGROUND_H
#define EXTRACTBACKGROUND_H
#include <armadillo>
#include <QFrame>
#include <QTimer>
#include "common.h"
namespace Ui {
class ExtractBackground;
}

class ExtractBackground : public QFrame
{
    Q_OBJECT
public:
    explicit ExtractBackground(
            std::vector<QWidget*>&inputs,
            QTimer& gl_timer,
            std::vector<arma::uvec>&labels,
            QWidget *parent = 0
            );
    bool configure(Config::Ptr);
    ~ExtractBackground();
signals:
    void message(QString,int);
protected:
    uint64_t current_frame_;
protected slots:
    void reset_extract_planes(){current_frame_=0;}
    void start_extract_planes();
    void finish_extract_planes();

    void extract_points();
    void unextract_points();
protected:

private:
    Ui::ExtractBackground *ui;
    std::vector<QWidget*>&inputs_;
    QTimer& gl_timer_;
    std::vector<arma::uvec>&labels_;
};

#endif // EXTRACTBACKGROUND_H
