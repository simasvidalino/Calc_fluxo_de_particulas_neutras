#ifndef DATAMATRICES_H_INCLUDED
#define DATAMATRICES_H_INCLUDED

#include <iostream>
#include <string>
#include <fstream>         //trabalhar com arquivo
#include <sstream>         //converter string para double, usar ostringstream para o titulo de saída
#include <vector>          //para acessar funcoes como eraser e shrink_to_fit() no vector
#include "LegendrePolynomial.h" //Para usar polinomio de legendre
#include "GausLegendreQuadrature.h"    //dados de quadraturas de Gauss-Legendre
#include "VariablesUsed.h"

using namespace std;

void construir_dados(string arquivotxt, dados_entrada &valor)
{
    double j;
    string k;
    ifstream dados_txt;
    vector <double> vetor;
    dados_txt.open(arquivotxt.c_str());


    if(!dados_txt)
    {
        cout<<"\nFalha na abertura do arquivo...\n"<<arquivotxt<<std::endl;
    }
    else
    {

        while(dados_txt){  //lê todo o arquivo
            dados_txt >> k;

            if(k[0] == '/'){ //pega somente linhas com texto (primeira palavra)
                dados_txt.ignore(1000,'\n');
            }else{
                istringstream(k) >> j;
                vetor.push_back(j);
            }
        }
    }

    dados_txt.close(); //fecha o arquivo txt

    /********************************************************************************************************

                           Distribuição de dados em matrizes

*********************************************************************************************************/

    if (vetor.empty())
        throw std::invalid_argument("Verifique se há erro no endereço do arquivo ou se ele não tem o formato padrão.");

    vetor.erase(vetor.begin(),vetor.begin()+9);
    static int i=0;

    //Ordem da Quadratura GL
    valor.n = vetor[i];
    i++;

    //Ordem de parada
    valor.ordem_parada = vetor[i];
    i++;

    //Maximo Numero de Iteracoes
    valor.iteracao = vetor[i];
    i++;

    //Grupos de Energia
    valor.G = vetor[i];
    i++;

    //Grau da Anisotropia do Espalhamento (L)
    valor.L = vetor[i];
    i++;

    //NR e NZ
    valor.n_R = vetor[i]; i++;

    valor.n_Z = vetor[i]; i++;

    //Tamanho de cada Regiao
    valor.TAM = new double [valor.n_R];
    for (int j=0;j<valor.n_R;j++){
        valor.TAM[j]=vetor[i];
        i++;
    }

    //Nodos por Regiao
    valor.n_nodos = new short int [valor.n_R];
    for (int j=0;j<valor.n_R;j++){
        valor.n_nodos[j]=vetor[i];
        i++;
    }

    //Periodicidade
    valor.periodicidade = vetor[i];
    i++;

    //Mapeamento
    valor.Map_R = new short int [valor.n_R];
    for (int j=0;j<valor.n_R;j++){
        valor.Map_R[j] = vetor[i];
        i++;
    }

    vetor.erase(vetor.begin(),vetor.begin()+i);
    vetor.shrink_to_fit() ;
    i=0;

    /*********************************************************************************************************

                                    Matrizes

*********************************************************************************************************/

    valor.fonte_g = new double*[valor.G];  //Fonte Fisica
    valor.s_t = new long double*[valor.G];      //Sigma total
    valor.s_s = new long double***[valor.G];    //Sigma Espalhamento

    for (int j=0;j<valor.G;j++){
        valor.fonte_g[j] = new double [valor.n_R];
        valor.s_s[j] = new long double**[valor.G];
        valor.s_t[j] = new long double[valor.n_Z];

        for (int k=0;k<valor.G;k++){
            valor.s_s[j][k] = new long double *[valor.n_Z];

            for (int l=0;l<valor.n_Z;l++){
                valor.s_s[j][k][l] = new long double [valor.L+1];}
        }
    }
    /**********************************************************************************************************/

    //Sigma total e Sigma de espalhamento
    for (int h=0;h<valor.n_Z;h++){
        for (int j=0;j<valor.G;j++){
            valor.s_t[j][h] = vetor[i];
            i++;
        }
        for (int k=0;k<valor.L+1;k++){
            for (int m=0;m<valor.G;m++){
                for (int n=0;n<valor.G;n++){
                    valor.s_s[m][n][h][k] = vetor[i];
                    i++;
                }
            }}}

    //Tipo de Condicoes de Contorno (Esq. Dir) (1-Prescrita. 2-Reflexiva)
    valor.tipo_ce = vetor[i];i++;                valor.tipo_cd = vetor[i]; i++;

    //Valor da condicao de contorno prescrita (Esq. Dir)
    valor.cceg = new double [valor.G];
    valor.ccdg = new double [valor.G];

    for (int j=0;j<valor.G;j++){
        valor.cceg[j] = vetor[i];
        i++;
        valor.ccdg[j] = vetor[i];
        i++;
    }

    //Fonte Fisica
    for (int j=0;j<valor.G;j++){
        for (int k=0;k<valor.n_R;k++){
            valor.fonte_g[j][k]=vetor[i];
            i++;
        }
    }
    vetor.clear();
    vetor.shrink_to_fit() ;

    /*********************************************************************************************************

                                           Dados calculados

*********************************************************************************************************/
    valor.CONTX = new int[valor.n_R]; //nodos por regiao acumulados
    valor.PASSO = new long double [valor.n_R]; //modulo entre nodos
    valor.TAM_TOTAL = 0; //comprimento total de x
    valor.NODOSX=0; // numero total de nodos por regiao

    for (int j=0;j<valor.n_R;j++){
        valor.CONTX[j] = valor.NODOSX+valor.n_nodos[j];
        valor.PASSO[j] = valor.TAM[j]/valor.n_nodos[j];
        valor.NODOSX = valor.CONTX[j];
    }
    valor.w = new double [valor.n];
    valor.mi = new double [valor.n];

    //fluxo angular e fluxo escalar

    valor.FLUXO_ANGULAR = new long double**[valor.G];
    valor.smgi          = new long double**[valor.G];
    valor.FLUXO_ESCALAR = new long double*[valor.G];

    valor.FLUXO_ANGULAR_DIREITA  = new long double**[valor.G];
    valor.FLUXO_ANGULAR_ESQUERDA = new long double**[valor.G];

    for (int g = 0; g < valor.G; g++)
    {
        valor.FLUXO_ANGULAR[g] = new long double*[(valor.NODOSX) + 1];
        valor.smgi[g]          = new long double*[valor.NODOSX      ];
        valor.FLUXO_ESCALAR[g] = new long double[(valor.NODOSX)  + 1];

        valor.FLUXO_ANGULAR_DIREITA[g]  = new long double*[(valor.NODOSX) + 1];
        valor.FLUXO_ANGULAR_ESQUERDA[g] = new long double*[(valor.NODOSX) + 1];

        for (int o = 0; o <= valor.NODOSX; o++)
        {
            valor.FLUXO_ANGULAR[g][o] = new long double[valor.n];

            valor.FLUXO_ANGULAR_DIREITA[g][o]   = new long double[valor.n];
            valor.FLUXO_ANGULAR_ESQUERDA[g][o]  = new long double[valor.n];
        }

        for (int o = 0; o < valor.NODOSX; o++)
        {
            valor.smgi[g][o] = new long double[valor.n];
        }
    }


    //Matriz com os polinômios de Legendre
    legendre_set ( valor.n,valor.mi,valor.w);
    valor.Mat_Legendre = new double *[valor.n];
    double *legendre_n = new double [valor.L+1];

    for (int n=0;n<valor.n;n++){
        valor.Mat_Legendre[n] = new double [valor.L + 1];
    }

    setlocale(LC_ALL,"portuguese");

    for (int n=0;n<valor.n;n++){
        Legendre::Pn(valor.L,valor.mi[n], legendre_n);
        for (int l=0;l<valor.L+1;l++){
            valor.Mat_Legendre[n][l] = legendre_n[l];
        }
    }
    i=0;

    delete [] legendre_n;
}

#endif // DATAMATRICES_H_INCLUDED
