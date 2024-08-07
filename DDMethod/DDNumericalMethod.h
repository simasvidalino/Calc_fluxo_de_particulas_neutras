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

void runDD();

inline void DD(dados_entrada &valor){

    /**************************************** método DD***********************************************************************

O Algoritmo

1- Colocar as condições de contorno nos fluxos angulares ψ
2- fazer a varredura para a direita e em seguida para a esquerda, obtendo fluxos angulares ψ em cada direção.
3- Calcular fluxo escalar
4- Testar convergência: a diferença entre o fluxo escalar calculado na iteração anterior com a atual deve ser menor do que
a tolerânicia para o processo parar.
5- Recalcular fonte de espalhamento e repetir todo o processo.
6- Caso a número de iterações chegue no valor declarado, o processo para.

************************************************************************************************************************/

    double tolerancia;
    tolerancia = 1/pow(10, valor.ordem_parada);

    //Inicialização do fluxo angular
    for(int g = 0; g < valor.G; g++)
    {
        for(int node = 0; node <= valor.NODOSX; node++)
        {
            for(int n = 0; n < valor.n; n++)
            {
                valor.FLUXO_ANGULAR[g][node][n] = 0.0;
            }
        }
    }

    //Condição de contorno pela esquerda e pela direita
    for(int g = 0; g < valor.G; g++)
    {
        for(int n = 0; n < valor.n; n++)
        {
            valor.FLUXO_ANGULAR[g][0][n]            = valor.cceg[g];
            valor.FLUXO_ANGULAR[g][valor.NODOSX][n] = valor.ccdg[g];
        }
    }

    //inicialização de smgi e fluxo escalar
    for(int g = 0; g < valor.G; g++)
    {
        for(int indexNode = 0; indexNode <= valor.NODOSX; indexNode++)
        {
            valor.FLUXO_ESCALAR[g][indexNode] = 0.0;
        }

        for(int indexNode = 0; indexNode < valor.NODOSX; indexNode++)
        {
            for(int n = 0; n < valor.n; n++)
            {
                valor.smgi[g][indexNode][n] = 0.0;
            }
        }
    }

    int iteracao = 0;
    long double aux;
    clock_t t1;

    while(iteracao < valor.iteracao)
    {

        int ni = 0;
        long double fi    = 0.0;
        long double Q     = 0.0;
        long double a     = 0.0;
        long double b     = 0.0;
        long double num   = 0.0;
        long double den   = 1.0;
        long double passo = 0.5;

        if(valor.tipo_ce == 2)
        {
            for(int g=0;g<valor.G;g++){
                int a=0;
                for(int o=(valor.n/2);o<valor.n;o++){
                    valor.FLUXO_ANGULAR[g][0][o] = valor.FLUXO_ANGULAR[g][0][a] ;
                    a++;
                }}}

        //varredura para a direita mi>0
        for(int g=0;g<valor.G;g++){
            ni = 0;
            for(int r=0;r<valor.n_R;r++)
            {
                Q     = valor.fonte_g[g][r];
                a     = 0.5 * valor.s_t[g][valor.Map_R[r] - 1];
                passo = valor.PASSO[r];
                for(int n=0;n<(valor.n_nodos[r]);n++)
                {     //n muda conforme mudamos de regiao, ou seja, cada regiao tem uma quantidade de nodos
                    for(int o=(valor.n/2);o<valor.n;o++)
                    {  //a varredura pela direita contempla a metadade das direções mi
                        b = valor.mi[o]/passo;
                        num = ((b - a) * valor.FLUXO_ANGULAR[g][ni][o]) + Q + valor.smgi[g][ni][o];
                        den = b + a;
                        fi = num/den;
                        valor.FLUXO_ANGULAR[g][ni + 1][o] = fi;
                    }
                    ni++;
                }
            }
        }

        if(valor.tipo_cd==2)
        {
            for(int g=0;g<valor.G;g++){
                int a=0;
                for(int o=0;o<valor.n/2;o++){
                    valor.FLUXO_ANGULAR[g][valor.NODOSX][o] = valor.FLUXO_ANGULAR[g][valor.NODOSX][a] ;
                    a++;
                }}}

        //varredura para a esquerda mi<0
        for(int g=0;g<valor.G;g++){
            ni = valor.NODOSX;
            for(int r=(valor.n_R-1);r>=0;r--)
            {
                Q     = valor.fonte_g[g][r];
                a     = 0.5*valor.s_t[g][valor.Map_R[r]-1];
                passo = valor.PASSO[r];

                for(int n=valor.n_nodos[r]-1;n>=0;n--)
                {
                    for(int o=0;o<(valor.n/2);o++)
                    {
                        b   = -valor.mi[o]/passo;
                        num = ((b - a)*valor.FLUXO_ANGULAR[g][ni][o]) + Q + valor.smgi[g][ni - 1][o];
                        den = b + a;
                        fi  = num/den;
                        valor.FLUXO_ANGULAR[g][ni - 1][o] = fi;
                    }
                    ni--;
                }
            }
        }

        //fluxo escalar
        long double somatorio1;
        aux=-1;

        for(int g=0;g<valor.G;g++){
            for(int n=0;n<=valor.NODOSX;n++){
                somatorio1 = 0;
                for(int o=0;o < valor.n;o++)
                {
                    somatorio1 = somatorio1 + (valor.FLUXO_ANGULAR[g][n][o])*(valor.w[o]);
                }
                somatorio1 = somatorio1 * 0.5;

                if(somatorio1 != 0)
                {
                    valor.MODULO = fabs((valor.FLUXO_ESCALAR[g][n] - somatorio1)/somatorio1);
                }
                else
                {
                    valor.MODULO = fabs(valor.FLUXO_ESCALAR[g][n] - somatorio1);
                }

                aux = max(aux,valor.MODULO); //salva-se o maior erro entre os fluxos escalares
                valor.FLUXO_ESCALAR[g][n] = somatorio1;

                if (std::isnan(valor.FLUXO_ESCALAR[g][n] ))
                {
                    std::cerr<<"Scalar flux is nan"<<std::endl;
                }
            }
        }

        /*cout<<"\n\n\n\n"<<aux<<"\n\n\n\n";
for(int g=0;g<valor.G;g++){
         cout<<" Grupo "<<valor.G<<"\n\n";
    for(int n=0;n<=valor.NODOSX;n++){
    cout<<" nodo "<<n<< "fluxo "<<valor.FLUXO_ESCALAR[g][n]<<endl;
    }}*/


        if(aux<tolerancia){//se a diferença entre os fluxos da iteração atual e anterior respeitar uma tolerância, então o programa para.
            break;
        }

        //recalcular smgi
        long double soma2,soma3;

        for(int g=0;g<valor.G;g++){
            ni=0;
            for(int r=0;r<valor.n_R;r++){
                for(int nodo=0;nodo<valor.n_nodos[r];nodo++){
                    for(int m=0;m<valor.n;m++){
                        soma3 = 0;
                        for(int g1=0;g1<valor.G;g1++)
                        {
                            for(int n=0;n<valor.n;n++)
                            {
                                soma2 = 0;
                                for(int l = 0; l < valor.L+1;l++)
                                {
                                    a = (2 * l) + 1;
                                    b = valor.Mat_Legendre[m][l] * valor.Mat_Legendre[n][l];
                                    soma2 = soma2 + (a * b * (valor.s_s[g][g1][(valor.Map_R[r])-1][l]) ) ;
                                }

                                soma3 = soma3 + ((valor.FLUXO_ANGULAR[g1][ni+1][n] + valor.FLUXO_ANGULAR[g1][ni][n])*valor.w[n]*soma2);

                                if (std::isnan(soma3))
                                {
                                    std::cerr<<"The Number is nan in smgi calculation"<<std::endl;
                                    soma3 = 0.0;
                                }
                            }
                        }
                        valor.smgi[g][ni][m] = 0.25*soma3; //o 0.25 surgiu de 0.5 da expressão somatorio2 vezes 0.5 do somatorio3
                        if (std::isnan(valor.smgi[g][ni][m] ))
                        {
                            std::cerr<<"The smgi is nan"<<std::endl;
                        }
                    }
                    ni++;
                }}}

        iteracao++;
    }

    t1 = clock(); //medição do tempo de execução

    valor.iteracaoFinal = iteracao;
    valor.tempoFinalDeProcessamento = (float)t1/CLOCKS_PER_SEC;
}

class DDMethod
{

public:
    static DDMethod *getInstance();
    void destroyInstance();

    void run(dados_entrada &valor);

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
