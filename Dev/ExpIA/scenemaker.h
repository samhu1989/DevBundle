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
            QWidget *parent = 0
            );
    bool configure(Config::Ptr);
    ~SceneMaker();
public slots:
    void save_to_inputs();
private:
    Ui::SceneMaker *ui;
    MeshListViewerWidget* lst_view_;
    std::vector<WidgetPtr>& mesh_views_;
    MeshBundle<DefaultMesh>::PtrList& inputs_;
};

#endif // SCENEMAKER_H
