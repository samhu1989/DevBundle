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
            MeshBundle<DefaultMesh>::PtrList& inputs,
            std::vector<arma::uvec>& labels,
            QObject *parent = 0
            );
    bool configure(Config::Ptr);
public:
    static std::vector<arma::umat> base_segment_list_;
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
    MeshBundle<DefaultMesh>::PtrList& inputs_;
    std::vector<arma::uvec>& labels_;
    arma::uword base_segment_i_;
    arma::uword base_segment_N_;
    Segmentation::RobustRegionDetection<DefaultMesh> cuts_;
};
#endif // ROBUSTCUT_H
