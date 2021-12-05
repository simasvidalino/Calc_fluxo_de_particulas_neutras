#include <iostream>
#include <string>
#include <fstream>         //trabalhar com arquivo
#include <ctime>           //para medir o tempo de execu��o
#include <cmath>           //para usar a fun��o pow
#include <sstream>         //converter string para double, usar ostringstream para o titulo de sa�da
#include <vector>          //para acessar funcoes como eraser e shrink_to_fit() no vector
#include <iomanip>         //para mostrar maior numero de casas decimais na tela
#include "DataMatrices.h" //struct para guardar vari�veis
#include "LegendrePolynomial.h" //Para usar polinomio de legendre
#include "GausLegendreQuadrature.h"    //dados de quadraturas de Gauss-Legendre
#include "DataMatrices.h" //construir dados de entrada
#include "DDNumericalMethold.h" //m�todo Diamond Difference
#include "VariablesUsed.h"

using namespace std;

int main()
{
    construir_dados("Arquivo_de_entrada.txt", valor);
    DD();

    return 0;
}
