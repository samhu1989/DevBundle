#ifndef UNIFYLABELCOLORSIZETHREAD_H
#define UNIFYLABELCOLORSIZETHREAD_H
#include <QThread>
#include "common.h"
class UnifyLabelColorSizeThread:public QThread
{
    Q_OBJECT
public:
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator InputIterator;
    typedef std::vector<arma::uvec>::iterator OutputIterator;
    UnifyLabelColorSizeThread(
            std::vector<MeshBundle<DefaultMesh>::Ptr>& inputs,
            std::vector<arma::uvec>& outputs
            ):inputs_(inputs),labels_(outputs)
    {
        setObjectName("UnifyLabelColorSizeThread");
    }
    bool configure(Config::Ptr config);
signals:
    void message(QString,int);
protected:
    void run();
private:
    std::vector<MeshBundle<DefaultMesh>::Ptr>& inputs_;
    std::vector<arma::uvec>& labels_;
    Config::Ptr config_;
};

#endif // UNIFYLABELCOLORSIZETHREAD_H
