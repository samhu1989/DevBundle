#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <QThread>
#include <QImage>
#include <armadillo>
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
public slots:
    void load_img(void);
    void load_annotation(void);

    void start_editing();
    void finish_editing();
protected:
    QImage input_img_;
    arma::uvec annotation_;
private:
    Ui::MainWindow *ui;
    QThread* edit_thread_;
};

#endif // MAINWINDOW_H
