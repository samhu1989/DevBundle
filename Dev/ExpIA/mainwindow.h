#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <QThread>
#include <common.h>
#include "objectmodel.h"
#include "voxelgraph.h"
using namespace OpenMesh;
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    typedef QWidget* WidgetPtr;
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

signals:
    void close_views();

public slots:
    void showInMdi(QWidget* w, Qt::WindowFlags flag = 0);
    void closeInMdi(QWidget* w);
    void showBox(int,MeshBundle<DefaultMesh>::Ptr);

protected slots:
    void configure();

    void open_inputs();
    void open_inputs(QStringList&);
    bool open_mesh(DefaultMesh&,const std::string&);
    void view_inputs();
    void removeView();

    void save_labels();
    void load_labels();

    void save_objects();
    void load_objects();

    void save_supervoxels();
    void load_supervoxels();

    void start_editing();
    void finish_editing();

    void showLab();
    void viewObj();
private:
    Ui::MainWindow *ui;
    std::vector<WidgetPtr> mesh_views_;
    MeshBundle<DefaultMesh>::PtrList inputs_;
    std::vector<arma::uvec> labels_;
    std::vector<ObjModel::Ptr> objects_;

    IO::Options io_opt_;
    Config::Ptr config_;
    QThread* edit_thread_;
};

#endif // MAINWINDOW_H
