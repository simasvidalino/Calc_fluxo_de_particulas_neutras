#pragma once

#include <functional>
#include "VariablesUsed.h"
using namespace std;

class DDMethod
{

public:
    static DDMethod *getInstance();
    void destroyInstance();

    void runDDMethodWithOMP(dados_entrada &data);
    void runDDMethodWithTwoThreads(dados_entrada &data);
    void runDDMethodWithOneThread(dados_entrada &data);

    void saveTxtFile(dados_entrada &data);

protected:
    virtual void initVariables();

    //It is not the one in the middle of node.
    virtual long double solveScalarFluxWithAbsoluteCriteria();
    virtual long double solveScalarFluxWithRelativeCriteria();

    virtual void leftScan();
    virtual void rightScan();

    virtual void leftScanReflective();
    virtual void rightScanReflective();

    virtual void calculateScatteringSourceSmgi();

    virtual void writeLog();

private:
    static DDMethod* myPointer;

    dados_entrada* values;
};
