#ifndef REGIONGROWRGBTHREAD_H
#define REGIONGROWRGBTHREAD_H
#include <QObject>
#include "common.h"
class RegionGrowRGBThread: public QObject
{
    Q_OBJECT
public:
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator InputIterator;
    typedef std::vector<arma::uvec>::iterator OutputIterator;
    explicit RegionGrowRGBThread(
            std::vector<MeshBundle<DefaultMesh>::Ptr>& inputs,
            std::vector<arma::uvec>& outputs,
            QObject *parent = 0
            );
    bool configure(Config::Ptr config);
public slots:
    void process(void);
signals:
    void message(QString,int);
    void finished(void);
private:
    std::vector<MeshBundle<DefaultMesh>::Ptr>& inputs_;
    std::vector<arma::uvec>& labels_;
    Config::Ptr config_;
    int verbose_;
};

#endif // REGIONGROWRGBTHREAD_H
