#ifndef ANNOTATOR_H
#define ANNOTATOR_H
#include <QFrame>
#include <armadillo>
#include <deque>
#include <memory>
#include <ext/hash_map>
#include "common.h"
namespace Ui {
class Annotator;
}

class Annotation
{
public:
    Annotation(
        arma::uvec indices,
        arma::uword tL,
        MeshBundle<DefaultMesh>::Ptr mesh,
        arma::uvec& label
    ):indices_(indices),tLabel_(tL),m_(mesh),target_(label)
    {
       ;
    }
    void apply(void)
    {
        sLabel_ = target_(indices_);
        target_(indices_).fill(tLabel_);
    }
    void undo(void)
    {
        target_(indices_) = sLabel_;
    }
    inline MeshBundle<DefaultMesh>::Ptr mesh_ptr()const{return m_;}
    inline const arma::uvec& label()const{return target_;}
private:
    MeshBundle<DefaultMesh>::Ptr m_;
    arma::uvec& target_;
    const arma::uvec indices_;
    arma::uvec sLabel_;
    const arma::uword tLabel_;
};

class Annotator : public QFrame
{
    Q_OBJECT
public:
    explicit Annotator(
            std::vector<QWidget*> &inputs,
            QTimer &gl_timer,
            std::vector<arma::uvec> &labels,
            QWidget *parent = 0
            );
    ~Annotator();
    bool configure(Config::Ptr);
    void init(void);
signals:
    void message(QString,int);
public slots:
    void undo(void);
protected slots:
    void pickLabel(void);
    void applyToPatch(void);
    void applyToPoint(void);
    void updateLabel(int label);
    void update_color(MeshBundle<DefaultMesh>::Ptr m,const arma::uvec&);
protected:
    void undo_annotation();
private:
    Ui::Annotator *ui;

    __gnu_cxx::hash_map<uint32_t,int> c2l_;
    __gnu_cxx::hash_map<int,uint32_t> l2c_;

    int max_queue_size_;
    std::deque<std::shared_ptr<Annotation>> annotation_deque_;

    std::vector<QWidget*>&inputs_;
    QTimer& gl_timer_;
    std::vector<arma::uvec>&labels_;
};

#endif // ANNOTATOR_H
