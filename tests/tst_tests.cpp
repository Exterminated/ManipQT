#include <QtTest>

// add necessary includes here

class tests : public QObject
{
    Q_OBJECT

public:
    tests();
    ~tests();

private slots:
    void initTestCase();
    void cleanupTestCase();
    void test_case1();

};

tests::tests()
{

}

tests::~tests()
{

}

void tests::initTestCase()
{

}

void tests::cleanupTestCase()
{

}

void tests::test_case1()
{
    double denominator_1, numerator_1, denominator_2, numerator_2;
    double fi = 0.072, OA1 =750, xmk = 282.35, ymk = 1191.748, zmk = -100.0, DK= -40.0, OK = 755.0, l1 = 1400.0, l4 = -40;
    denominator_1 = (2.0*OA1*cos(fi)*(ymk+OA1*sin(fi))+2.0*OA1*sin(fi)*(zmk-OA1*cos(fi)))*(l1-sqrt(pow(ymk+OA1*sin(fi),2.0)+pow(zmk-OA1*cos(fi),2.0)));
    numerator_1 = sqrt(pow(ymk+OA1*sin(fi),2.0)+pow(zmk-OA1*cos(fi),2.0)+xmk*xmk);
    denominator_2 = (l4-sqrt(pow(DK-OA1*cos(fi),2.0)+pow(OK-OA1*sin(fi),2.0)))*(2*OA1*sin(fi)*(DK-OA1*cos(fi))-2.0*OA1*cos(fi)*(OK-OA1*sin(fi)));
    numerator_2 = sqrt(pow(DK-OA1*cos(fi),2.0)+pow(OK-OA1*sin(fi),2.0));

}

QTEST_APPLESS_MAIN(tests)

#include "tst_tests.moc"
