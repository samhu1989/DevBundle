#ifndef SCENEMAKER_H
#define SCENEMAKER_H
#include <QWidget>
#include "configure.h"
#include "MeshListViewerWidget.h"
namespace Ui {
class SceneMaker;
}

class SceneMaker : public QWidget
{
    Q_OBJECT
public:
    typedef QWidget* WidgetPtr;
    explicit SceneMaker(
            std::vector<WidgetPtr>& views,
            MeshBundle<DefaultMesh>::PtrList& inputs,
            std::vector<arma::uvec>& labels,
            QWidget *parent = 0
            );
    bool configure(Config::Ptr);
    ~SceneMaker();
signals:
    void message(QString,int);
    void view_input(WidgetPtr);
public slots:
    void save_to_inputs();
    void save_to_input(DefaultMesh&);
    void save_to_label(arma::uvec&);
    void save_rt();
    void reset_rt();
    void update_rt(arma::fmat R,arma::fvec t,int idx);
private:
    Ui::SceneMaker *ui;
    MeshListViewerWidget* lst_view_;
    std::vector<WidgetPtr>& mesh_views_;
    MeshBundle<DefaultMesh>::PtrList& inputs_;
    std::vector<arma::uvec>& labels_;
    std::vector<arma::fmat> R_lst_;
    std::vector<arma::fvec> t_lst_;
};

#endif // SCENEMAKER_H
