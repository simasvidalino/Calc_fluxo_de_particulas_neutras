#ifndef ESTRUTURAS_H_INCLUDED
#define ESTRUTURAS_H_INCLUDED

struct dados_entrada{

/* Dados fornecidos pelo usuário
****************************************************************************************
* n            -> ordem da quadratura
* w            -> peso quadratura
* mi           -> direção (raiz polinomio de Legendre)
* L            -> Grau de espalhamento
* G            -> Número de grupos de energia
* Z            -> Numero de zonas materiais
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
double **fonte_g,**c_c;
long double ****s_s;
short int n,L,G,Z,n_R,n_Z,tipo_ce,tipo_cd, ordem_parada;
double periodicidade;
short int *n_nodos,*Map_R;
int iteracao;


/* Dados calculados
*******************************************************************************************
* CONTX        -> nodos por região acumulativos, para ter referência de localização X
* ERRO         -> critério de parada calculado a partir da Ordem_parada
* PASSO        -> dimensão do nodo espacial por região na direção X
*TAM_TOTAL     -> comprimento total do domínio X
*TOTAL_NODOS   -> total de nodos no domínio X
*FLUXO_ANGULAR1-> Fluxo de neutrons em cada nodo, grupo de energia e direção a ser calculado (anterior)
*FLUXO_ANGULAR2-> Fluxo de neutrons em cada nodo, grupo de energia e direção a ser calculado (posterior)
*******************************************************************************************/

long double MODULO;
double TAM_TOTAL;
long double *PASSO;
double*legendre_n;
int *CONTX;
long double ***FLUXO_ANGULAR;
long double **FLUXO_ESCALAR;
long double***smgi;
int NODOSX;
double **Mat_Legendre;

}valor;


#endif // ESTRUTURAS_H_INCLUDED
