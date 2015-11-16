#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <QThread>
#include <common.h>
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

protected slots:
    void configure();

    void open_inputs();
    void open_inputs(QStringList&);
    bool open_mesh(DefaultMesh&,const std::string&);
    void view_inputs();
    void removeView();

    void save_labels();
    void load_labels();

    void start_editing();
    void finish_editing();
private:
    Ui::MainWindow *ui;
    std::vector<MeshBundle<DefaultMesh>::Ptr> inputs_;
    std::vector<WidgetPtr> mesh_views_;
    std::vector<arma::uvec> labels_;
    IO::Options io_opt_;
    Config::Ptr config_;
    QThread* edit_thread_;
};

#endif // MAINWINDOW_H
