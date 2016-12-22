#ifndef JRCS_H
#define JRCS_H
#include "common.h"
#include "objectmodel.h"
#include <QObject>
#include <memory>
#include <armadillo>
#include "jrcsinitbase.h"
#include "jrcsview.h"
class JRCSWork : public QObject
{
    Q_OBJECT
public:
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> MeshList;
    typedef std::vector<arma::uvec> LabelList;
    typedef std::vector<ObjModel::Ptr> ModelList;
    typedef std::shared_ptr<arma::mat> DMatPtr;
    typedef std::vector<DMatPtr> DMatPtrLst;
    explicit JRCSWork(
            MeshList& inputs,
            LabelList& labels,
            ModelList& objects,
            QObject *parent = 0
            );
    ~JRCSWork(){}
signals:
    void end();
    void message(QString,int);
public:
    bool configure(Config::Ptr config);
public slots:
    void Init_SI_HKS();
    void debug_SI_HKS();
    void Init_Bernolli();
protected:
    void Init_Bernolli_a(arma::uvec&);
    void Init_Bernolli_b(arma::uvec&);
public:
    static void set_opt_aoni(JRCSView* w);
    static void set_opt_aopt(JRCSView* w);
    static void set_opt_spectrum(JRCSView* w);
    static void set_opt_bilateral(JRCSView* w);
    static bool init_optimize(JRCSView* w);
private:
    MeshList& inputs_;
    LabelList& labels_;
    ModelList& objects_;
    Config::Ptr config_;
    static DMatPtrLst alpha_ptrlst_;
    static arma::fvec obj_prob_;
};

#endif // JRCS_H
