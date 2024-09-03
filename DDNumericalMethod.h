#pragma once

#include <iostream>
#include <string>
#include <fstream>         //trabalhar com arquivo
#include <ctime>           //para medir o tempo de execução
#include <cmath>           //para usar a função pow
#include <sstream>         //converter string para double, usar ostringstream para o titulo de saída
#include <vector>          //para acessar funcoes como eraser e shrink_to_fit() no vector
#include <iomanip>         //para mostrar maior numero de casas decimais na tela
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
