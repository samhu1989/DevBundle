#ifndef MAKEINCOMPLETE_H
#define MAKEINCOMPLETE_H
#include "common.h"
#include <QObject>
#include <QString>
#include <QList>
#include <random>
class MakeIncomplete : public QObject
{
    Q_OBJECT
public:
    explicit MakeIncomplete(MeshBundle<DefaultMesh>::PtrList inputs,
            std::vector<arma::uvec> &labels,
            QList<float> ratio,
            QString path,
            QObject *parent = 0
            );
    bool configure(Config::Ptr);
signals:
    void end(void);
    void message(QString,int);
public slots:
    void process(void);
protected:
    void saveByRatio(MeshBundle<DefaultMesh>::Ptr in, arma::uvec &li);
private:
    MeshBundle<DefaultMesh>::PtrList inputs_;
    std::vector<arma::uvec>& labels_;
    QList<float> ratio_;
    QString path_;
    std::default_random_engine gen_;
};

#endif // MAKEINCOMPLETE_H
