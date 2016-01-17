#ifndef UNIFYLABELMANNUAL_H
#define UNIFYLABELMANNUAL_H
#include <armadillo>
#include <QFrame>
#include "common.h"
#include <QKeyEvent>
#include <QPalette>
#include <QMdiSubWindow>
#include "patchpairview.h"
#include "MeshViewerWidget.h"
#include <QTimer>
namespace Ui {
class UnifyLabelMannual;
}

class UnifyLabelMannual : public QFrame
{
    Q_OBJECT
public:
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> InputType;
    typedef std::vector<arma::uvec> OutputType;
    typedef std::shared_ptr<DefaultMesh> MeshPtr;
    explicit UnifyLabelMannual(InputType&input,OutputType&output,QWidget*parent=0);
    ~UnifyLabelMannual();
    bool configure(Config::Ptr);
    void initLater(){timer.start(20);}
signals:
    void message(QString,int);
public slots:
    void init(void);
protected:
    void keyPressEvent(QKeyEvent*);

    void frameNext();
    void frameLast();

    void movePatchNext();
    void movePatchLast();
    void patchNext();
    void patchDelete();

    void showHelp();

    void reloadFrame();
    void showPatches();
    void updateLabel();

    void updateObjects();
private:
    Ui::UnifyLabelMannual *ui;
    InputType& inputs_;
    OutputType& labels_;
    std::vector<PatchPairView*> pair_views_;
    std::vector<MeshPtr> objects_;
    std::vector<MeshPtr> patches_;
    arma::uvec patch_label_;
    size_t current_frame_;
    size_t current_patch_;
    Config::Ptr config_;
    QTimer timer;
};

#endif // UNIFYLABELMANNUAL_H
