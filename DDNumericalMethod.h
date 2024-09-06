#pragma once

//#include "LegendrePolynomial.h" //Para usar polinomio de legendre
//#include "GausLegendreQuadrature.h"    //dados de quadraturas de Gauss-Legendre
#include "VariablesUsed.h"
using namespace std;

class DDMethod
{

public:
    static DDMethod *getInstance();
    void destroyInstance();

    void runDDMethodWithTwoThreads(dados_entrada &valor);
    void runDDMethodWithOneThread(dados_entrada &valor);

    void saveTxtFile(dados_entrada &valor);

protected:
    virtual void leftScan();
    virtual void rightScan();
    virtual void initVariables();
    virtual void updateAngularFluxMatrix();
    virtual void updateAngularFluxMatrixLeft();
    virtual void updateAngularFluxMatrixRight();

private:
    static DDMethod* myPointer;

    dados_entrada* values;
};
