#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <QThread>
#include <common.h>
#include "objectmodel.h"
#include "voxelgraph.h"
#include <QTimer>
#include "tests.h"
#include <QDir>
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
    void object_updated();
    void keyPressSignal(QKeyEvent*);

public slots:
    void showInMdi(QWidget* w, Qt::WindowFlags flag = 0);
    void closeInMdi(QWidget* w);
    void showBox(int,MeshBundle<DefaultMesh>::Ptr);
    void notify_object_updated(){
        std::cerr<<"notifying"<<std::endl;
        emit object_updated();
    }

public slots:
    void save_labels(QString dirName=QString());
    void load_labels();

    void save_base_segs(QString dirName=QString());
    void load_base_segs();

    void save_objects(QString dirName=QString());
    void load_objects();

    void save_supervoxels(QString dirName=QString());
    void load_supervoxels();

    void save_cluster(QString dirName=QString());
    void load_cluster();

    void save_vox_index_picked(QString dirName=QString());
    void load_vox_index_picked();

    void save_pts_index_picked(QString dirName=QString());
    void load_pts_index_picked();

    void save_Pix_Order_Functor(QString dirName=QString());
    void save_Vox_Order_Functor(QString dirName=QString());

    void save_cube_color(QString path=QString());
    void load_cube_color(QString path=QString());

    void view_input(WidgetPtr widget);
    void calculate_iou(QString dir0=QString(),QString dir1=QString());
    void calculate_fit(QString dir=QString());

    void save_XYZRGBL_MAT(QString dirName=QString());
protected slots:
    void configure();

    void open_inputs();
    void open_inputs(QDir& dir);
    void open_inputs(QStringList&);
    void save_aligned();
    bool open_mesh(DefaultMesh&,const std::string&);
    void view_inputs();

    void removeView();

    void save_scenes();
    void save_object_layout(const std::string&);

    void save_scene_layout(const std::string&);
    void save_scene_model(const std::string&);

    void keyPressEvent(QKeyEvent* event){emit keyPressSignal(event);}

    void make_scene();
    void start_editing();
    void update_object();
    void finish_editing();
    void finish_editing_ui();

    void remove_zero_label();

    void showLab();
    void viewObj();
    void showSVColor();
    void showFeature();
    void showIndex();
    void showSpectralFunc();

    void goOver();

    void custom_color_from_cube();

protected:
    QAction* getActionByText(const QString& txt);
    void save_XYZRGBL_MAT(MeshBundle<DefaultMesh>::Ptr,arma::uvec&,QString);

protected slots:
    void LAPACKE_dggsvd_test(void){TEST::LAPACKE_dggsvd_test();}
    void Inside_BBox_test(void){TEST::Inside_BBox_test();}
    void agd_test(void){TEST::agd_test();}
    void jrcs_plate_test(void){TEST::jrcs_plate_test();}
    void jrcs_cube_test(void){TEST::jrcs_cube_test();}

private:
    Ui::MainWindow *ui;
    std::vector<WidgetPtr> mesh_views_;
    MeshBundle<DefaultMesh>::PtrList inputs_;
    std::vector<arma::uvec> labels_;
    std::vector<ObjModel::Ptr> objects_;
    arma::mat feature_base_;
    arma::mat feature_centers_;

    IO::Options io_opt_;
    Config::Ptr config_;
    QThread* edit_thread_;
    QWidget* edit_widget_;
    QTimer gl_timer;
};

#endif // MAINWINDOW_H
