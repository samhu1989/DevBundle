#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <QThread>
#include <common.h>
#include "objectmodel.h"
#include "voxelgraph.h"
#include <QTimer>
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
    void save_aligned();
    bool open_mesh(DefaultMesh&,const std::string&);
    void view_inputs();
    void removeView();

    void save_labels();
    void load_labels();

    void save_objects();
    void load_objects();

    void save_supervoxels();
    void load_supervoxels();

    void save_scenes();
    void save_object_layout(const std::string&);
    void save_scene_layout(const std::string&);
    void save_scene_model(const std::string&);

    void start_editing();
    void finish_editing();

    void remove_zero_label();

    void showLab();
    void viewObj();
    void showSVColor();
private:
    Ui::MainWindow *ui;
    std::vector<WidgetPtr> mesh_views_;
    MeshBundle<DefaultMesh>::PtrList inputs_;
    std::vector<arma::uvec> labels_;
    std::vector<ObjModel::Ptr> objects_;

    IO::Options io_opt_;
    Config::Ptr config_;
    QThread* edit_thread_;
    QTimer gl_timer;
};

#endif // MAINWINDOW_H
