#ifndef VARIABLESUSED_H_INCLUDED
#define VARIABLESUSED_H_INCLUDED

struct dados_entrada{

/* Dados fornecidos pelo usu�rio
****************************************************************************************
* n            -> ordem da quadratura
* w            -> peso quadratura
* mi           -> dire��o (raiz polinomio de Legendre)
* L            -> Grau de espalhamento
* G            -> N�mero de grupos de energia
* Z            -> Numero de zonas materiais
* n_R          -> Numero de regi�es na dire��o x
* n_Z          -> Numero de zonas materiais
* p            -> Periodicidade DD-Periodicidade para mostrar resultados
* TAM          -> Comprimento de cada regiao
* n_nodos      -> N�mero de nodos em cada regi�o na dire��o x
* Map_R        -> Mapeamento das regi�oes
* s_s          -> Se��o de choque de espalhamento (L,G e Z)
* s_t          -> Se��o de choque total (G e R)
* c_c          -> Condi��o de contorno
* iteracao     -> N�mero m�ximo de itera��es
* tipo_ce      -> Tipo da condi��o de contorno da esquerda
* tipo_cd      -> Tipo da condi��o de contorno da direita
* fonte_g      -> Intensidade da fonte por regi�o no problema f�sico (1D) Multigrupo
* Ordem_parada -> Ordem do valor a atingir para a parada. Ex: Ordem_parada=7, ent�o o erro � 10^7
* cceg - ccdg  -> valor do fluxo f�sico incidente para o problema f�sico Cond Cont Prescrita por Grupo
* periodicidade-> valor, em nodos/cm, que os dados estar�o no arquivo de s�ida

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
* CONTX        -> nodos por regi�o acumulativos, para ter refer�ncia de localiza��o X
* ERRO         -> crit�rio de parada calculado a partir da Ordem_parada
* PASSO        -> dimens�o do nodo espacial por regi�o na dire��o X
*TAM_TOTAL     -> comprimento total do dom�nio X
*TOTAL_NODOS   -> total de nodos no dom�nio X
*FLUXO_ANGULAR1-> Fluxo de neutrons em cada nodo, grupo de energia e dire��o a ser calculado (anterior)
*FLUXO_ANGULAR2-> Fluxo de neutrons em cada nodo, grupo de energia e dire��o a ser calculado (posterior)
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


#endif // VARIABLESUSED_H_INCLUDED
