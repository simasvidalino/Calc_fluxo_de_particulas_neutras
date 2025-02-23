#pragma once

#include <iostream>
#include <memory>
#include <vector>


struct dados_entrada
{

    /* Dados fornecidos pelo usuário
****************************************************************************************
* n            -> ordem da quadratura
* w            -> peso quadratura
* mi           -> direção (raiz polinomio de Legendre)
* L            -> Grau de espalhamento
* G            -> Número de grupos de energia
* n_R          -> Numero de regiões na direção x
* n_Z          -> Numero de zonas materiais
* p            -> Periodicidade DD-Periodicidade para mostrar resultados
* TAM          -> Comprimento de cada regiao
* n_nodos      -> Número de nodos em cada região na direção x
* Map_R        -> Mapeamento das regiãoes
* s_s          -> Seção de choque de espalhamento (L,G e Z)
* s_t          -> Seção de choque total (G e R)
* c_c          -> Condição de contorno
* iteracao     -> Número máximo de iterações
* tipo_ce      -> Tipo da condição de contorno da esquerda
* tipo_cd      -> Tipo da condição de contorno da direita
* fonte_g      -> Intensidade da fonte por região no problema físico (1D) Multigrupo
* Ordem_parada -> Ordem do valor a atingir para a parada. Ex: Ordem_parada=7, então o erro é 10^7
* cceg - ccdg  -> valor do fluxo físico incidente para o problema físico Cond Cont Prescrita por Grupo
* periodicidade-> valor, em nodos/cm, que os dados estarão no arquivo de sáida

****************************************************************************************/
    double *w, *mi,*TAM,*cceg,*ccdg;
    long double **s_t;
    double **fonte_g;
    long double ****s_s;
    short int n = 0, L = 0, G = 0, n_R = 0, n_Z = 0, tipo_ce = 0, tipo_cd = 0, ordem_parada = 0;
    double periodicidade = 0.0;
    short int *n_nodos,*Map_R;
    int iteracao = 0;
    int iteracaoFinal = 0;
    float tempoFinalDeProcessamento = 0.0;

    /* Dados calculados
*******************************************************************************************
* CONTX        -> nodos por região acumulativos, para ter referência de localização X
* ERRO         -> critério de parada calculado a partir da Ordem_parada
* PASSO        -> dimensão do nodo espacial por região na direção X
*TAM_TOTAL     -> comprimento total do domínio X
*TOTAL_NODOS   -> total de nodos no domínio X
*FLUXO_ANGULAR1-> Fluxo de neutrons em cada nodo, grupo de energia e direção a ser calculado (anterior)
*FLUXO_ANGULAR2-> Fluxo de neutrons em cada nodo, grupo de energia e direção a ser calculado (posterior)
*s_ab          -> A seção de choque de absorção em cada zona representa a probabilidade de um nêutron ser
*                 absorvido pelo material. A seção de choque de absorção para um determinado grupo de
*                 energia pode ser calculada subtraindo a seção de choque de espalhamento da seção de choque total
*                 desse grupo.
*******************************************************************************************/
    long double MODULO;
    double TAM_TOTAL;
    long double *PASSO;
    int *CONTX;
    long double ***FLUXO_ANGULAR;
    long double ***FLUXO_ANGULAR_DIREITA;
    long double ***FLUXO_ANGULAR_ESQUERDA;

    long double **FLUXO_ESCALAR;
    long double***smgi;
    int NODOSX;
    double **Mat_Legendre;

    ~dados_entrada()
    {

        std::cout<<"Call delete matrices"<<std::flush;

        auto safeDelete = [](auto& ptr) {
            if (ptr != nullptr) {
                delete[] ptr;
                ptr = nullptr;
            }
        };

        safeDelete(w);
        safeDelete(mi);
        safeDelete(TAM);
        safeDelete(cceg);
        safeDelete(ccdg);
        safeDelete(n_nodos);
        safeDelete(Map_R);
        safeDelete(PASSO);
        safeDelete(CONTX);

        // Desalocar memória para s_s
        if (s_s != nullptr)
        {
            for (int j = 0; j < G; ++j) {
                for (int k = 0; k < G; ++k) {
                    for (int l = 0; l < n_Z; ++l) {
                        safeDelete(s_s[j][k][l]);
                    }
                    safeDelete(s_s[j][k]);
                }
                safeDelete(s_s[j]);
            }
            safeDelete(s_s);
        }


        // Desalocar memória para s_t
        if (s_t != nullptr)
        {
            for (int j = 0; j < G; ++j)
            {
                safeDelete(s_t[j]);
            }
            safeDelete(s_t);
        }

        // Desalocar memória para fonte_g
        if (fonte_g != nullptr) {
            for (int j = 0; j < G; ++j)
            {
                safeDelete(fonte_g[j]);
            }
            safeDelete(fonte_g);
        }

        // Desalocar memória para Condição de Contorno
        if (fonte_g != nullptr)
        {
            for (int j = 0; j < G; ++j)
            {
                safeDelete(fonte_g[j]);
            }
            safeDelete(fonte_g);
        }

        // Desalocar memória para FLUXO_ANGULAR
        if (FLUXO_ANGULAR != nullptr)
        {
            for (int g = 0; g < G; ++g) {
                for (int o = 0; o <= NODOSX; ++o)
                {
                    safeDelete(FLUXO_ANGULAR[g][o]);
                    safeDelete(FLUXO_ANGULAR_DIREITA[g][o]);
                    safeDelete(FLUXO_ANGULAR_ESQUERDA[g][o]);
                }
                safeDelete(FLUXO_ANGULAR[g]);
                safeDelete(FLUXO_ANGULAR_ESQUERDA[g]);
                safeDelete(FLUXO_ANGULAR_DIREITA[g]);
            }

            safeDelete(FLUXO_ANGULAR);
            safeDelete(FLUXO_ANGULAR_DIREITA);
            safeDelete(FLUXO_ANGULAR_ESQUERDA);
        }

        if (FLUXO_ANGULAR_DIREITA != nullptr)
        {
            for (int g = 0; g < G; ++g) {
                for (int o = 0; o <= NODOSX; ++o)
                {

                    safeDelete(FLUXO_ANGULAR_DIREITA);
                }
                safeDelete(FLUXO_ANGULAR_DIREITA[g]);
            }

            safeDelete(FLUXO_ANGULAR_DIREITA);
        }

        if (FLUXO_ANGULAR_ESQUERDA != nullptr) {
            for (int g = 0; g < G; ++g) {
                for (int o = 0; o <= NODOSX; ++o) {
                    safeDelete(FLUXO_ANGULAR_ESQUERDA);
                }
                safeDelete(FLUXO_ANGULAR_ESQUERDA[g]);
            }
            safeDelete(FLUXO_ANGULAR_ESQUERDA);
        }

        // Desalocar memória para FLUXO_ESCALAR
        if (FLUXO_ESCALAR != nullptr) {
            for (int g = 0; g < G; ++g) {
                safeDelete(FLUXO_ESCALAR[g]);
            }
            safeDelete(FLUXO_ESCALAR);
        }

        // Desalocar memória para smgi
        if (smgi != nullptr) {
            for (int g = 0; g < G; ++g) {
                for (int o = 0; o < NODOSX; ++o) {
                    safeDelete(smgi[g][o]);
                }
                safeDelete(smgi[g]);
            }
            safeDelete(smgi);
        }

        // Desalocar memória para Mat_Legendre
        if (Mat_Legendre != nullptr) {
            for (int index = 0; index < n; ++index) {
                safeDelete(Mat_Legendre[index]);
            }
            safeDelete(Mat_Legendre);
        }

        //17 variables
    }
};







