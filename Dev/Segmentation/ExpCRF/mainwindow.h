#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <QThread>
#include <QImage>
#include <armadillo>
#include <QTimer>
#include "common.h"
#include "MeshLabelViewerWidget.h"
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
protected slots:
    void showInMdi(QWidget* w, Qt::WindowFlags flag = 0);

    void load_img(void);
    void load_annotation(void);

    void load_mesh(void);
    void read_mesh(const QString& filename);
    void view_mesh(void);

    void start_editing();
    void finish_editing();

protected:
    bool open_mesh(DefaultMesh&,const std::string&);

protected:
    QImage input_img_;
    std::shared_ptr<arma::uvec> annotation_;
    QHash<arma::uword,arma::uword> colortolabel_;

    MeshBundle<DefaultMesh>::Ptr input_mesh_;
    std::shared_ptr<MeshLabelViewerWidget> input_mesh_view_;
    QTimer gl_timer_;
private:
    Ui::MainWindow *ui;
    QThread* edit_thread_;
};

#endif // MAINWINDOW_H
