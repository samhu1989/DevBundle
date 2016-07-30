#ifndef SORT_AGD_H
#define SORT_AGD_H
#include <QObject>
#include "common.h"
class Sort_AGD:public QObject
{
    Q_OBJECT
public:
    typedef MeshBundle<DefaultMesh>::PtrList InputList;
    Sort_AGD(
            MeshBundle<DefaultMesh>::PtrList& inputs,
            QObject* parent=0
            );
    bool configure(Config::Ptr);
public slots:
    void process(void);
signals:
    void finished();
    void message(QString,int);
protected:
    void sort(const arma::vec& agd,MeshBundle<DefaultMesh>& m);
private:
    MeshBundle<DefaultMesh>::PtrList& inputs_;
};

#endif // SORT_AGD_H
