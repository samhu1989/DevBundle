#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <QMessageBox>
#include <QMenuBar>
#include <QGLWidget>
#include <QThread>
#include <QKeyEvent>
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
    void keyPressEvent(QKeyEvent*);
    void start_registration(void);
    void finish_registration(void);
private:
    Ui::MainWindow *ui;
    QThread* alg_thread;
};



#endif // MAINWINDOW_H
