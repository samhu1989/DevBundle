#include <QString>
#include <QtTest>
#include <armadillo>
#include "iocore.h"
class Io_Test : public QObject
{
    Q_OBJECT
public:
    Io_Test();
private Q_SLOTS:
    void writemat();
    void writefmat();
    void wirtecmat();
};

Io_Test::Io_Test()
{
}

void Io_Test::writemat()
{
    arma::mat a = {{43.5, 5432.434, 3.32},{0.0001, 88834.0, 0.0}};
    MATIO::save_to_matlab(a,"./debug/a.mat");

    arma::vec b = {43.5, 5432.434, 3.32, 0.0001, 88834.0};
    MATIO::save_to_matlab(b,"./debug/b.mat");

    arma::rowvec c = {43.5, 5432.434, 3.32, 0.0001, 88834.0};
    MATIO::save_to_matlab(c,"./debug/c.mat");
}

void Io_Test::writefmat()
{
    arma::fmat fa = {{43.5, 5432.434, 3.32},{0.0001, 88834.0, 0.0}};
    MATIO::save_to_matlab(fa,"./debug/fa.mat");

    arma::fvec fb = {43.5, 5432.434, 3.32, 0.0001, 88834.0};
    MATIO::save_to_matlab(fb,"./debug/fb.mat");

    arma::frowvec fc = {43.5, 5432.434, 3.32, 0.0001, 88834.0};
    MATIO::save_to_matlab(fc,"./debug/fc.mat");
}

void Io_Test::wirtecmat()
{
    arma::fmat ca = {{255,255,255},{255,0,255},{0,0,0}};
    MATIO::save_to_matlab(ca,"./debug/ca.mat");
}

QTEST_APPLESS_MAIN(Io_Test)
#include "tst_io_test.moc"
