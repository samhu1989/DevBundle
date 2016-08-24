#ifndef NCUT_H
#define NCUT_H
#include "common.h"
#include "segmentationcore.h"
#include <QObject>
class NCut : public QObject
{
    Q_OBJECT
public:
    explicit NCut(
            MeshBundle<DefaultMesh>::PtrList& inputs,
            std::vector<arma::uvec>& labels,
            QObject *parent = 0
            );
signals:
    void end();
    void message(QString,int);
public slots:
    bool configure(Config::Ptr);
    void process();
    void debug_convexity();
    void debug_W();
private:
    MeshBundle<DefaultMesh>::PtrList& inputs_;
    std::vector<arma::uvec>& labels_;
    Segmentation::NormalizedCuts<DefaultMesh> cuts_;
};

#endif // NCUT_H
