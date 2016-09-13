#ifndef ROBUSTCUT_H
#define ROBUSTCUT_H
#include "common.h"
#include "segmentationcore.h"
#include <QObject>
#include <QKeyEvent>
class RobustCut : public QObject
{
    Q_OBJECT
public:
    explicit RobustCut(
            QImage&,
            arma::uvec& labels,
            QObject *parent = 0
            );
    bool configure(Config::Ptr);
public:
    static arma::umat base_segment_;
signals:
    void end();
    void message(QString,int);
public slots:
    void base_segments();
    void keyPressEvent(QKeyEvent*);
    void consensus_segment();
    void show_base_segment();
protected:
    void next_base_segment();
    void last_base_segment();
    void exit_base_segment();
private:
    const QImage& inputs_;
    arma::uvec& labels_;
    arma::uword base_segment_i_;
    arma::uword base_segment_N_;
    Segmentation::RobustRegionDetection<DefaultMesh> cuts_;
};
#endif // ROBUSTCUT_H
