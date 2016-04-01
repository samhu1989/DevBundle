#include <QString>
#include <QtTest>
#include "optimizationcore.h"
#include "lbfgs.h"
class OptimizationTest : public QObject
{
    Q_OBJECT

public:
    OptimizationTest();
private Q_SLOTS:
    void LBFGS_Case1();
};

OptimizationTest::OptimizationTest()
{

}

class ExampleEnergy: public Optimization::EnergyFunction {
public:
    virtual arma::vec initialValue() {
        return arma::vec( 2 , arma::fill::zeros );
    }
    virtual double gradient( const arma::vec & x, arma::vec & dx ) {
        double fx = (x[0] - 1)*(x[0] - 6) + (x[0] - 4)*(x[1] - 2)*(x[0] - 4)*(x[1] - 2) + x[1]*x[1];
        dx[0] = 2*x[0] - 7 + 2*(x[1] - 2)*(x[0] - 4)*(x[1] - 2);
        dx[1] = 2*x[1] + 2*(x[0] - 4)*(x[1] - 2)*(x[0] - 4);
        return fx;
    }
};

void OptimizationTest::LBFGS_Case1()
{
    ExampleEnergy e;
    Optimization::LBFGS opt;
    arma::vec m = opt.minimize(e,0,false);
    arma::vec g(m.n_rows);
    std::cout<<"Result:"<<m<<std::endl;
    std::cout<<"Minimized :"<<e.gradient(m,g)<<std::endl;
    QVERIFY2(true, "Failure");
}

QTEST_APPLESS_MAIN(OptimizationTest)

#include "tst_optimizationtest.moc"
