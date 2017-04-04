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
    explicit SceneMaker(QWidget *parent = 0);
    bool configure(Config::Ptr);
    ~SceneMaker();
private:
    Ui::SceneMaker *ui;
    MeshListViewerWidget* lst_view_;
};

#endif // SCENEMAKER_H
