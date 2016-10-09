#include <QString>
#include <QtTest>

class Io_testTest : public QObject
{
    Q_OBJECT

public:
    Io_testTest();

private Q_SLOTS:
    void testCase1();
};

Io_testTest::Io_testTest()
{
}

void Io_testTest::testCase1()
{
    QVERIFY2(true, "Failure");
}

QTEST_APPLESS_MAIN(Io_testTest)

#include "tst_io_testtest.moc"
