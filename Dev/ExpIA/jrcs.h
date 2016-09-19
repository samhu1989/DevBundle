#ifndef JRCS_H
#define JRCS_H
#include <QObject>
#include <memory>
#include <armadillo>
class JRCS : public QObject
{
    Q_OBJECT
public:
    typedef std::shared_ptr<arma::fmat> MatPtr;
    typedef std::vector<MatPtr> MatPtrLst;
    explicit JRCS(QObject *parent = 0);
signals:
    void end();
public slots:
    void Init_SI_HSK();
public:
    static void optimization();
private:
    static MatPtrLst alpha_ptrlst_;
};

#endif // JRCS_H
