#ifndef LABELCOMPACTOR_H
#define LABELCOMPACTOR_H
#include <QObject>
#include <armadillo>
#include "common.h"
class LabelCompactor : public QObject
{
    Q_OBJECT
public:
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> MeshList;
    typedef std::vector<arma::uvec> LabelList;
    explicit LabelCompactor(
            MeshList& inputs,
            LabelList& labels,
            QObject *parent = 0
            );
    bool configure(Config::Ptr);
signals:
    void end();
    void message(QString,int);
public slots:
    void process();
private:
    MeshList& inputs_;
    LabelList& labels_;
};

#endif // LABELCOMPACTOR_H
