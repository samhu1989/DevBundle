#include <QString>
#include <QtTest>
#include "optimizationcore.h"
#include "lbfgs.h"
#include "sdp.h"
class OptimizationTest : public QObject
{
    Q_OBJECT

public:
    OptimizationTest();
private Q_SLOTS:
    void LBFGS_Case1();
    void SDP_Example();
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

class Example_SDP: public Optimization::SDP
{
public:
    Example_SDP()
    {
        std::vector<arma::mat> C;
        C.reserve(3);
        arma::mat c0 = {{2,1},{1,2}};
        C.push_back(c0);
        arma::mat c1 = {{3,0,1},{0,2,0},{1,0,3}};
        C.push_back(c1);
        C.emplace_back(2,1,arma::fill::zeros);//input diag block
        setC(C);
        std::vector<std::vector<arma::mat>> As;
        As.reserve(2);
        As.emplace_back(3);
        arma::mat a00 = {{3,1},{1,3}};
        As[0][0] = a00;
        As[0][1].clear();
        arma::vec a02 = {1,0};
        As[0][2] = a02;
        As.emplace_back(3);
        As[1][0] = arma::vec(2,arma::fill::zeros);//input diag block
        arma::mat a11 = {{3,0,1},{0,4,0},{1,0,5}};
        As[1][1] = a11;
        arma::vec a12 = {0,1};
        As[1][2] = a12;
        setAs(As);
        ey_<<0.75<<1.00<<arma::endr;
        arma::vec value;
        value<<0.25<<-0.25<<-0.25<<0.25<<2.00<<2.00<<0.75<<1.00<<arma::endr;
        arma::umat location;
        location<<0<<1<<0<<1<<3<<4<<5<<6<<arma::endr
                <<0<<0<<1<<1<<3<<4<<5<<6<<arma::endr;
        eZ_ = arma::sp_mat(location,value);
        value.clear();
        value<<0.125<<0.125<<0.125<<0.125<<0.667<<arma::endr;
        location.clear();
        location<<0<<1<<0<<1<<2<<arma::endr
                <<0<<0<<1<<1<<2<<arma::endr;
        eX_ = arma::sp_mat(location,value);
    }
    bool test()
    {
        if(!init())return false;
        if(!solve())return false;
        arma::vec y;
        arma::sp_mat Z;
        arma::sp_mat X;
        gety(y);
        getZ(Z);
        getX(X);
        bool pass = true;
        if(!arma::approx_equal(y,ey_,"absdiff",1e-10))
        {
            pass = false;
            std::cerr<<"result y:"<<std::endl;
            std::cerr<<y<<std::endl;
            std::cerr<<"expect y:"<<std::endl;
            std::cerr<<ey_<<std::endl;
        }
        if(!arma::approx_equal(Z,eZ_,"absdiff",1e-10))
        {
            pass = false;
            std::cerr<<"result Z:"<<std::endl;
            std::cerr<<Z<<std::endl;
            std::cerr<<"expect Z:"<<std::endl;
            std::cerr<<eZ_<<std::endl;
        }
        if(!arma::approx_equal(X,eX_,"absdiff",1e-10))
        {
            pass = false;
            std::cerr<<"result X:"<<std::endl;
            std::cerr<<X<<std::endl;
            std::cerr<<"expect X:"<<std::endl;
            std::cerr<<eX_<<std::endl;
        }
        return pass;
    }
private:
    arma::vec ey_;
    arma::sp_mat eX_;
    arma::sp_mat eZ_;
};

void OptimizationTest::SDP_Example()
{
    Example_SDP sdp;
    QVERIFY2(sdp.test(),"SDP Failed");
}

QTEST_APPLESS_MAIN(OptimizationTest)

#include "tst_optimizationtest.moc"
