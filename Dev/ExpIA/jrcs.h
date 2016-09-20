#ifndef JRCS_H
#define JRCS_H
#include "common.h"
#include "objectmodel.h"
#include <QObject>
#include <memory>
#include <armadillo>
class JRCSWork : public QObject
{
    Q_OBJECT
public:
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> MeshList;
    typedef std::vector<arma::uvec> LabelList;
    typedef std::vector<ObjModel::Ptr> ModelList;
    typedef std::shared_ptr<arma::fmat> MatPtr;
    typedef std::vector<MatPtr> MatPtrLst;
    explicit JRCSWork(
            MeshList& inputs,
            LabelList& labels,
            ModelList& objects,
            QObject *parent = 0
            );
signals:
    void end();
    void message(QString,int);
public:
    bool configure(Config::Ptr config);
public slots:
    void Init_SI_HSK();
public:
    static void optimization();
private:
    MeshList& inputs_;
    LabelList& labels_;
    ModelList& objects_;
    static MatPtrLst alpha_ptrlst_;
};

#endif // JRCS_H
