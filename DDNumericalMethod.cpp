#include "DDNumericalMethod.h"

#include <thread>
#include <mutex>
#include <iostream>
#include <string>
#include <fstream>         //trabalhar com arquivo
#include <ctime>           //para medir o tempo de execução
#include <cmath>           //para usar a função pow
#include <sstream>         //converter string para double, usar ostringstream para o titulo de saída
#include <iomanip>         //para mostrar maior numero de casas decimais na tela


std::mutex mtxl;
std::mutex mtxr;

DDMethod* DDMethod::myPointer = nullptr;

DDMethod *DDMethod::getInstance()
{
    if (nullptr == myPointer)
        myPointer = new DDMethod;

    return myPointer;
}

void DDMethod::destroyInstance()
{
    if (nullptr != myPointer)
    {
        delete myPointer;
        myPointer = nullptr;
    }
}

void DDMethod::runDDMethodWithTwoThreads(dados_entrada &valor)
{
    values = &valor;

    initVariables();

    int ni = 0;
    long double fi    = 0.0;
    long double a     = 0.0;
    long double b     = 0.0;
    long double num   = 0.0;
    int iteration = 0;
    long double aux;

    const auto start_time = std::chrono::high_resolution_clock::now();

    double tol;
    tol = 1.0/pow(10, values->ordem_parada);

    while (iteration < values->iteracao)
    {
        std::thread righthread(&DDMethod::rightScan, this);
        std::thread leftThread(&DDMethod::leftScan, this);

        leftThread.join();
        righthread.join();

        //fluxo escalar
        long double somatorio1;
        aux=-1;

        for(int g=0;g<values->G;g++)
        {
            for(int n=0;n<=values->NODOSX;n++)
            {
                somatorio1 = 0;
                for(int o=0;o < values->n;o++)
                {
                    somatorio1 = somatorio1 + (values->FLUXO_ANGULAR[g][n][o])*(values->w[o]);
                }
                somatorio1 = somatorio1 * 0.5;

                if(somatorio1 != 0)
                {
                    values->MODULO = fabs((values->FLUXO_ESCALAR[g][n] - somatorio1)/somatorio1);
                }
                else
                {
                    values->MODULO = fabs(values->FLUXO_ESCALAR[g][n] - somatorio1);
                }

                aux = max(aux,values->MODULO); //salva-se o maior erro entre os fluxos escalares
                values->FLUXO_ESCALAR[g][n] = somatorio1;

                if (std::isnan(values->FLUXO_ESCALAR[g][n] ))
                {
                    std::cerr<<"Scalar flux is nan"<<std::endl;
                }
            }
        }

        // cout<<"\n\n\n\n"<<aux<<"\n\n\n\n";
        // for(int g=0;g<values->G;g++){
        //     cout<<" Grupo "<<values->G<<"\n\n";
        //     for(int n=0;n<=values->NODOSX;n++){
        //         cout<<" nodo "<<n<< "fluxo "<<values->FLUXO_ESCALAR[g][n]<<endl;
        //     }}


        if(aux<tol)
        {//se a diferença entre os fluxos da iteração atual e anterior respeitar uma tolerância, então o programa para.
            break;
        }

        //recalcular smgi
        long double soma2,soma3;

        for(int g=0;g<values->G;g++){
            ni=0;
            for(int r=0;r<values->n_R;r++){
                for(int nodo=0;nodo<values->n_nodos[r];nodo++){
                    for(int m=0;m<values->n;m++){
                        soma3 = 0;
                        for(int g1=0;g1<values->G;g1++)
                        {
                            for(int n=0;n<values->n;n++)
                            {
                                soma2 = 0;
                                for(int l = 0; l < values->L+1;l++)
                                {
                                    a = (2 * l) + 1;
                                    b = values->Mat_Legendre[m][l] * values->Mat_Legendre[n][l];
                                    soma2 = soma2 + (a * b * (values->s_s[g][g1][(values->Map_R[r])-1][l]) ) ;
                                }

                                soma3 = soma3 + ((values->FLUXO_ANGULAR[g1][ni+1][n] + values->FLUXO_ANGULAR[g1][ni][n])*values->w[n]*soma2);

                                if (std::isnan(soma3))
                                {
                                    std::cerr<<"The Number is nan in smgi calculation"<<std::endl;
                                    soma3 = 0.0;
                                }
                            }
                        }
                        values->smgi[g][ni][m] = 0.25*soma3; //o 0.25 surgiu de 0.5 da expressão somatorio2 vezes 0.5 do somatorio3
                        if (std::isnan(values->smgi[g][ni][m] ))
                        {
                            std::cerr<<"The smgi is nan"<<std::endl;
                        }
                    }
                    ni++;
                }}}

        iteration++;
    }

    const auto end_time               = std::chrono::high_resolution_clock::now();
    values->tempoFinalDeProcessamento = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
                                        / 1000.0f;
    values->iteracaoFinal = iteration;
}

void DDMethod::runDDMethodWithOneThread(dados_entrada &valor)
{
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
    values = &valor;

    long double tolerancia;
    tolerancia = 1 / pow(10, valor.ordem_parada);

    initVariables();

    int iteracao = 0;
    long double aux;
    const auto start_time = std::chrono::high_resolution_clock::now();

    // for (int i = 0; i < valor.n; ++i)
    // {
    //     std::cout << std::fixed << std::setprecision(15) << "i = " << i << ", mi_hardcoded = " << valor.mi[i]
    //               << ", w_hardcoded = " << valor.w[i] << std::endl;
    // }

    while (iteracao < valor.iteracao.load(std::memory_order_relaxed))
    {
        int ni                   = 0;

        long double fi           = 0.0;
        long double Q            = 0.0;
        long double a            = 0.0;
        long double b            = 0.0;
        long double num          = 0.0;
        long double den          = 1.0;
        long double passo        = 0.5;

        short int stRegionIndice = 0;

        if (valor.tipo_ce == 2)
        {
            for (int g = 0; g < valor.G; g++)
            {
                int a = 0;
                for (int o = (valor.n / 2); o < valor.n; o++)
                {
                    valor.FLUXO_ANGULAR[g][0][o] = valor.FLUXO_ANGULAR[g][0][a];
                    a++;
                }
            }
        }

        //varredura para a direita mi > 0
        // psi(i+1) = [(b - a) * psi(i) + Q + smgi] / (b + a)
        for (int g = 0; g < valor.G; g++)
        {
            // ni representa o índice global do nodo atual na malha 1D. Neste caso, ele é crescente
            ni = 0;
            for (int r = 0; r < valor.n_R; r++)
            {
                stRegionIndice = valor.Map_R[r] - 1;

                Q              = valor.fonte_g[g][r];
                a              = 0.5 * valor.s_t[g][stRegionIndice];
                passo          = valor.PASSO[r];

                for (int n = 0; n < (valor.n_nodos[r]); n++)
                {
                    //n muda conforme mudamos de regiao, ou seja, cada regiao tem uma quantidade de nodos
                    for (int o = (valor.n / 2); o < valor.n; o++)
                    {
                        //a varredura pela direita contempla a metadade das direções mi
                        b   = valor.mi[o] / passo;
                        num = ((b - a) * valor.FLUXO_ANGULAR[g][ni][o]) + Q + valor.smgi[g][ni][o];
                        den = b + a;
                        fi  = num / den;
                        valor.FLUXO_ANGULAR[g][ni + 1][o] = fi;
                    }
                    ni++;
                }
            }
        }

        if (valor.tipo_cd == 2)
        {
            for (int g = 0; g < valor.G; g++)
            {
                int a = 0;
                for (int o = 0; o < valor.n / 2; o++)
                {
                    valor.FLUXO_ANGULAR[g][valor.NODOSX][o] = valor.FLUXO_ANGULAR[g][valor.NODOSX][a];
                    a++;
                }
            }
        }

        //varredura para a esquerda mi < 0
        //psi(i-1) = [(b - a) * psi(i) + Q + smgi] / (b + a)
        for (int g = 0; g < valor.G; g++)
        {
            // ni representa o índice global do nodo atual na malha 1D. Neste caso, ele é decrescente
            ni = valor.NODOSX;
            for (int r = (valor.n_R - 1); r >= 0; r--)
            {
                Q     = valor.fonte_g[g][r];
                a     = 0.5 * valor.s_t[g][valor.Map_R[r] - 1];
                passo = valor.PASSO[r];

                for (int n = valor.n_nodos[r] - 1; n >= 0; n--)
                {
                    for (int o = 0; o < (valor.n / 2); o++)
                    {
                        b   = -valor.mi[o] / passo;
                        num = ((b - a) * valor.FLUXO_ANGULAR[g][ni][o]) + Q + valor.smgi[g][ni - 1][o];
                        den = b + a;
                        fi  = num / den;
                        valor.FLUXO_ANGULAR[g][ni - 1][o] = fi;
                    }
                    ni--;
                }
            }
        }

        //fluxo escalar nos nodos (não é o médio)
        long double somatorio1;
        aux = -1;

        for (int g = 0; g < valor.G; g++)
        {
            for (int n = 0; n <= valor.NODOSX; n++)
            {
                somatorio1 = 0;
                for (int o = 0; o < valor.n; o++)
                {
                    somatorio1 = somatorio1 + (valor.FLUXO_ANGULAR[g][n][o]) * (valor.w[o]);
                }
                somatorio1 = somatorio1 * 0.5;

                if (somatorio1 != 0)
                {
                    valor.MODULO = fabs((valor.FLUXO_ESCALAR[g][n] - somatorio1)/* / somatorio1*/);
                }
                else
                {
                    valor.MODULO = fabs(valor.FLUXO_ESCALAR[g][n] - somatorio1);
                }

                aux                       = max(aux, valor.MODULO); //salva-se o maior erro entre os fluxos escalares
                valor.FLUXO_ESCALAR[g][n] = somatorio1;

                if (std::isnan(valor.FLUXO_ESCALAR[g][n]))
                {
                    std::cerr << "Scalar flux is nan" << std::endl;
                }
            }
        }

        // cout<<"\n\n\n\n"<<aux<<"\n\n\n\n";
        // for(int g=0;g<valor.G;g++){
        //     cout<<" Grupo "<<valor.G<<"\n\n";
        //     for(int n=0;n<=valor.NODOSX;n++){
        //         cout<<" nodo "<<n<< "fluxo "<<valor.FLUXO_ESCALAR[g][n]<<endl;
        //     }}

        // std::ofstream logFile;

        // logFile.open("logfile.txt", std::ios::app); // Abre o arquivo em modo de adição

        // if (!logFile.is_open())
        // {
        //     std::cerr << "Erro ao abrir o arquivo de log!" << std::endl;
        //     return;
        // }

        // logFile << "Inicializando variáveis...\n";

        // logFile << "\n" << aux << "<" << tolerancia << "\n\n\n\n";

        // for (int g = 0; g < valor.G; g++)
        // {
        //     logFile << " Grupo " << g << "\n\n";
        //     for (int n = 0; n <= valor.NODOSX; n++)
        //     {
        //         logFile << " nodo " << n << " fluxo " << valor.FLUXO_ESCALAR[g][n] << std::endl;
        //     }
        // }

        if (aux < tolerancia)
        {
            //se a diferença entre os fluxos da iteração atual e anterior respeitar uma tolerância, então o programa para.
            break;
        }

        //recalcular smgi
        long double soma2 = 0;
        long double soma3 = 0;

        std::vector<std::vector<std::vector<long double>>>
            legendreProj(valor.n, std::vector<std::vector<long double>>(valor.n, std::vector<long double>(valor.L + 1)));

        // Precalcula a parte (2l+1) * P_l(mu_m) * P_l(mu_n) para todas as direcoes e ordens
        for (int m = 0; m < valor.n; m++)
        {
            for (int n = 0; n < valor.n; n++)
            {
                for (int l = 0; l <= valor.L; l++)
                {
                    legendreProj[m][n][l] = (2 * l + 1) * valor.Mat_Legendre[m][l] * valor.Mat_Legendre[n][l];
                }
            }
        }

        for (int g = 0; g < valor.G; g++)
        {
            // ni representa o índice global do nodo atual na malha 1D.
            int ni = 0;

            for (int r = 0; r < valor.n_R; r++)
            {
                int z = valor.Map_R[r] - 1;

                for (int nodo = 0; nodo < valor.n_nodos[r]; nodo++)
                {
                    for (int m = 0; m < valor.n; m++)
                    {
                        long double soma3 = 0.0;

                        for (int g1 = 0; g1 < valor.G; g1++)
                        {
                            for (int n = 0; n < valor.n; n++)
                            {
                                //Fluxo interpolado deveria ser a média, mas para termos de otimização, a divisão por dois está abaixo
                                long double fluxoInterp = valor.FLUXO_ANGULAR[g1][ni][n] + valor.FLUXO_ANGULAR[g1][ni + 1][n];
                                long double pesoDir     = valor.w[n];
                                long double termoEspalhamentoAnisotropico = 0.0;

                                for (int l = 0; l <= valor.L; l++)
                                {
                                    // componenteAngularEspalhamento = soma_l( (2l+1) * P_l(mu_m) * P_l(mu_n) * sigma_s_l(g1 ->g) )
                                    termoEspalhamentoAnisotropico += legendreProj[m][n][l] * valor.s_s[g][g1][z][l];
                                }

                                soma3 += fluxoInterp * pesoDir * termoEspalhamentoAnisotropico;
                            }
                        }

                        //0.5 de Legendre vezes 0.5 da média do fluxo angular (fluxoInterp/2)
                        valor.smgi[g][ni][m] = 0.25 * soma3;

                        if (std::isnan(valor.smgi[g][ni][m]))
                        {
                            std::cerr << "The smgi is nan" << std::endl;
                        }
                    }

                    ni++;
                }
            }
        }

        iteracao++;
    }

    const auto end_time               = std::chrono::high_resolution_clock::now();
    values->tempoFinalDeProcessamento = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
                                        / 1000.0f;
    values->iteracaoFinal           = iteracao;

}

void DDMethod::saveTxtFile(dados_entrada &valor)
{
    ofstream saida_dados;
    string titulo;
    ostringstream auxiliar;

    // Gera o nome do arquivo baseado nos parâmetros
    auxiliar << "Saida_Dados_R" << valor.n_R << "_G" << valor.G << "_L" << valor.L << "_N" << valor.n << "_Nod" << valor.NODOSX << ".txt";
    titulo = auxiliar.str();

    // Abre o arquivo para escrita
    saida_dados.open(titulo);

    // Cabeçalho do arquivo
    saida_dados << string(160, '*') << endl;
    saida_dados << setw(90) << "Saída de dados" << endl;
    saida_dados << string(160, '*') << endl;

    // Informações gerais sobre o processamento
    saida_dados << "Número de iter:\t\t" << valor.iteracaoFinal << endl;
    saida_dados << "Tempo(s):\t\t" << valor.tempoFinalDeProcessamento / CLOCKS_PER_SEC << endl; // Converte o tempo para segundos
    saida_dados << string(160, '_') << endl;
    saida_dados << string(160, '*') << endl;

    // Título da seção de fluxo angular final
    saida_dados << "Fluxo angular final\t\titer: " << valor.iteracaoFinal << endl;
    saida_dados << string(160, '*') << endl;

    // Cabeçalhos da tabela de resultados
    saida_dados << "Grupo\t\t";
    valor.TAM_TOTAL = 0;

    // Calcula o comprimento total e imprime a primeira linha de distâncias
    for (int r = 0; r <= valor.n_R; r++) {
        valor.TAM_TOTAL += valor.TAM[r];
    }

    double t = 0;
    while (t <= valor.TAM_TOTAL) {
        saida_dados << setw(15)  << t << "cm";
        t += valor.periodicidade;
    }
    saida_dados << endl;

    saida_dados << string(160, '*') << endl;

    // Escrita dos fluxos escalares para cada grupo
    for (int g = 0; g < valor.G; g++) {
        saida_dados << g + 1 << "\t\t"; // Coluna do grupo
        int nod = 0;
        int steps = static_cast<int>(valor.TAM_TOTAL / valor.periodicidade);

        for (int i = 0; i <= steps; ++i)
        {
            saida_dados << setw(15) << scientific << setprecision(6) << valor.FLUXO_ESCALAR[g][nod];
            nod += (valor.NODOSX * valor.periodicidade) / valor.TAM_TOTAL;
        }
        saida_dados << endl;
    }

    saida_dados.close();
}

void DDMethod::leftScan()
{
    int ni = 0;
    long double fi    = 0.0;
    long double Q     = 0.0;
    long double a     = 0.0;
    long double b     = 0.0;
    long double num   = 0.0;
    long double den   = 1.0;
    long double passo = 0.5;

    std::unique_lock<std::mutex> lock(mtxr);

    if(values->tipo_cd == 2)
    {
        for(int g=0;g<values->G;g++){
            int a=0;
            for(int o=0;o<values->n/2;o++){
                values->FLUXO_ANGULAR_ESQUERDA[g][values->NODOSX][o] =
                    values->FLUXO_ANGULAR[g][values->NODOSX][a] ;
                a++;
            }}}

    lock.unlock();

    //varredura para a esquerda mi<0
    for(int g=0;g < values->G;g++)
    {
        ni = values->NODOSX;
        for(int r=(values->n_R-1);r>=0;r--)
        {
            Q     = values->fonte_g[g][r];
            a     = 0.5*values->s_t[g][values->Map_R[r]-1];
            passo = values->PASSO[r];

            for(int n=values->n_nodos[r]-1;n>=0;n--)
            {
                for(int o=0;o<(values->n/2);o++)
                {
                    b   = -values->mi[o]/passo;
                    num = ((b - a)*values->FLUXO_ANGULAR_ESQUERDA[g][ni][o]) + Q + values->smgi[g][ni - 1][o];
                    den = b + a;
                    fi  = num/den;
                    values->FLUXO_ANGULAR_ESQUERDA[g][ni - 1][o] = fi;
                }
                ni--;
            }
        }
    }

    updateAngularFluxMatrixLeft();
}

void DDMethod::rightScan()
{
    int ni = 0;
    long double fi    = 0.0;
    long double Q     = 0.0;
    long double a     = 0.0;
    long double b     = 0.0;
    long double num   = 0.0;
    long double den   = 1.0;
    long double passo = 0.5;

    std::unique_lock<std::mutex> lock(mtxr);

    if(values->tipo_ce == 2)
    {
        for(int g=0;g<values->G;g++){
            int a=0;
            for(int o=(values->n/2);o<values->n;o++){
                values->FLUXO_ANGULAR_DIREITA[g][0][o] = values->FLUXO_ANGULAR[g][0][a] ;
                a++;
            }}}

    lock.unlock();

    //varredura para a direita mi>0
    for(int g=0;g<values->G;g++)
    {
        ni = 0;
        for(int r=0;r<values->n_R;r++)
        {
            Q     = values->fonte_g[g][r];
            a     = 0.5 * values->s_t[g][values->Map_R[r] - 1];
            passo = values->PASSO[r];

            for(int n=0;n<(values->n_nodos[r]);n++)
            {     //n muda conforme mudamos de regiao, ou seja, cada regiao tem uma quantidade de nodos
                for(int o=(values->n/2);o<values->n;o++)
                {  //a varredura pela direita contempla a metadade das direções mi
                    b = values->mi[o]/passo;
                    num = ((b - a) * values->FLUXO_ANGULAR_DIREITA[g][ni][o]) + Q + values->smgi[g][ni][o];
                    den = b + a;
                    fi = num/den;
                    values->FLUXO_ANGULAR_DIREITA[g][ni + 1][o] = fi;
                }
                ni++;
            }
        }
    }

    updateAngularFluxMatrixRight();
}

void DDMethod::initVariables()
{
    //Inicialização do fluxo angular
    for(int g = 0; g < values->G; g++)
    {
        for(int node = 0; node <= values->NODOSX; node++)
        {
            for(int n = 0; n < values->n; n++)
            {
                values->FLUXO_ANGULAR[g][node][n] = 0.0;

                values->FLUXO_ANGULAR_DIREITA[g][node][n]  = 0.0;
                values->FLUXO_ANGULAR_ESQUERDA[g][node][n] = 0.0;
            }
        }
    }

    //Condição de contorno pela esquerda e pela direita
    for(int g = 0; g < values->G; g++)
    {
        for(int n = 0; n < values->n; n++)
        {
            values->FLUXO_ANGULAR[g][0][n]              = values->cceg[g];
            values->FLUXO_ANGULAR[g][values->NODOSX][n] = values->ccdg[g];

            values->FLUXO_ANGULAR_DIREITA[g][0][n]               = values->cceg[g];
            values->FLUXO_ANGULAR_ESQUERDA[g][values->NODOSX][n] = values->ccdg[g];
        }
    }

    //inicialização de smgi e fluxo escalar
    for(int g = 0; g < values->G; g++)
    {
        for(int indexNode = 0; indexNode <= values->NODOSX; indexNode++)
        {
            values->FLUXO_ESCALAR[g][indexNode] = 0.0;
        }

        for(int indexNode = 0; indexNode < values->NODOSX; indexNode++)
        {
            for(int n = 0; n < values->n; n++)
            {
                values->smgi[g][indexNode][n] = 0.0;
            }
        }
    }
}

void DDMethod::updateAngularFluxMatrix()
{
    for (int g = 0; g < values->G; g++) {
        for (int nodeIndex = 0; nodeIndex <= values->NODOSX; nodeIndex++)
        {
            // Atualização das direções de fluxo para a direita (mi > 0)
            for (int o = 0; o < values->n / 2; o++) {  // Apenas metade das direções (mi > 0)
                values->FLUXO_ANGULAR[g][nodeIndex][o] = values->FLUXO_ANGULAR_DIREITA[g][nodeIndex][o];
            }

            // Atualização das direções de fluxo para a esquerda (mi < 0)
            for (int o = values->n / 2; o < values->n; o++)
            {  // A outra metade das direções (mi < 0)
                values->FLUXO_ANGULAR[g][nodeIndex][o] = values->FLUXO_ANGULAR_ESQUERDA[g][nodeIndex][o - (values->n / 2)];
            }
        }
    }
}

void DDMethod::updateAngularFluxMatrixLeft()
{
    for (int g = 0; g < values->G; g++)
    {
        for (int nodeIndex = 0; nodeIndex <= values->NODOSX; nodeIndex++)
        {
            // Atualização das direções de fluxo para a esquerda (mi < 0)
            for (int o = values->n / 2; o < values->n; o++)
            {  // A outra metade das direções (mi < 0)
                values->FLUXO_ANGULAR[g][nodeIndex][o] = values->FLUXO_ANGULAR_ESQUERDA[g][nodeIndex][o - (values->n / 2)];
            }
        }
    }
}

void DDMethod::updateAngularFluxMatrixRight()
{
    for (int g = 0; g < values->G; g++)
    {
        for (int nodeIndex = 0; nodeIndex <= values->NODOSX; nodeIndex++)
        {
            // Atualização das direções de fluxo para a direita (mi > 0)
            for (int o = 0; o < values->n / 2; o++) {  // Apenas metade das direções (mi > 0)
                values->FLUXO_ANGULAR[g][nodeIndex][o] = values->FLUXO_ANGULAR_DIREITA[g][nodeIndex][o];
            }
        }
    }
}

