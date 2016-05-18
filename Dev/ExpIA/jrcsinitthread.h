#ifndef JRCSINITTHREAD_H
#define JRCSINITTHREAD_H
#include "common.h"
#include <QObject>
#include "hierarchicalization.h"
class JRCSInitThread:public QObject
{
    Q_OBJECT
signals:
    void finished();
    void message(QString,int);
    void showbox(int,MeshBundle<DefaultMesh>::Ptr);
public:
    JRCSInitThread(
            MeshBundle<DefaultMesh>::PtrList& inputs,
            std::vector<arma::uvec>& labels,
            QObject* parent = 0
            );
    bool configure(Config::Ptr);
public slots:
    void process();
protected:
private:
    Hierarchicalization seg_;
    MeshBundle<DefaultMesh>::PtrList& inputs_;
    std::vector<arma::uvec>& labels_;
};

#endif // JRCSINITTHREAD_H
