#pragma once

#include <functional>
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
