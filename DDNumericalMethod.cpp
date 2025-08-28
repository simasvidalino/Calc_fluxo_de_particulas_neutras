#include "DDNumericalMethod.h"

#include <thread>
#include <mutex>
#include <iostream>
#include <string>
#include <fstream>         //trabalhar com arquivo
#include <cmath>           //para usar a função pow
#include <sstream>         //converter string para double, usar ostringstream para o titulo de saída
#include <iomanip>         //para mostrar maior numero de casas decimais na tela

#include <atomic>
#include <omp.h>
#include <future>

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

void DDMethod::runDDMethodWithOMP(dados_entrada &data)
{
#ifdef USE_OPENMP
    const int max_threads = omp_get_max_threads();
    std::cout << "OpenMP enabled, max_threads = " << max_threads<< " threads\n"<<std::endl;

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
    values = &data;
    long double somatorio1;
    long double aux = -1;

    initVariables();
    long double iteration_max_error = 0.0;

    const long double tolerancia = 1.0 / std::pow(10.0L, values->ordem_parada);
    const auto &w = values->w;
    const auto &n = values->n;

    // Pré-computar funções
    auto rightScanFunc = (values->tipo_ce == 2) ? &DDMethod::rightScanReflective : &DDMethod::rightScan;
    auto leftScanFunc = (values->tipo_cd == 2) ? &DDMethod::leftScanReflective : &DDMethod::leftScan;

    int iteracao = 0;
    long double max_error = 0.0;
    const auto start_time = std::chrono::high_resolution_clock::now();

    // Criar pool de threads uma única vez
#pragma omp parallel
    {
        while (iteracao < values->iteracao.load(std::memory_order_relaxed))
        {
#pragma omp sections
            {
#pragma omp section
                (this->*rightScanFunc)();

#pragma omp section
                (this->*leftScanFunc)();
            }
             #pragma omp for collapse(2)
            // Restante do cálc

            // Barreira implícita - espera tasks completarem
            //fluxo escalar nos nodos (não é o médio)

            for (int g = 0; g < values->G; g++)
            {
                for (int n = 0; n <= values->NODOSX; n++)
                {
                    somatorio1 = 0;

                    for (int o = 0; o < values->n; o++)
                    {
                        somatorio1 = somatorio1 + (values->FLUXO_ANGULAR[g][n][o]) * (values->w[o]);
                    }

                    somatorio1 = somatorio1 * 0.5;

                    if (somatorio1 != 0)
                    {
                        values->MODULO = fabs((values->FLUXO_ESCALAR[g][n] - somatorio1) / somatorio1);
                    }
                    else
                    {
                        values->MODULO = fabs(values->FLUXO_ESCALAR[g][n] - somatorio1);
                    }

                    aux = max(aux, values->MODULO); //salva-se o maior erro entre os fluxos escalares

                    if (   ( somatorio1 < 0.0)
                        || ( std::isnan( somatorio1 ) ) )
                    {
                        //std::cerr << "Scalar flux is nan or negative " << somatorio1 << std::endl;
                        //values->FLUXO_ESCALAR[g][n] = 0.0;
                    }
                    else
                    {
                        //std::cerr << "Scalar flux is Positive " << somatorio1 << std::endl;

                        values->FLUXO_ESCALAR[g][n] = somatorio1;
                    }
                }
            }

            // cout<<"\n\n\n\n"<<aux<<"\n\n\n\n";
            // for(int g=0;g<values->G;g++){
            //     cout<<" Grupo "<<values->G<<"\n\n";
            //     for(int n=0;n<=values->NODOSX;n++){
            //         cout<<" nodo "<<n<< "fluxo "<<values->FLUXO_ESCALAR[g][n]<<endl;
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

            // for (int g = 0; g < values->G; g++)
            // {
            //     logFile << " Grupo " << g << "\n\n";
            //     for (int n = 0; n <= values->NODOSX; n++)
            //     {
            //         logFile << " nodo " << n << " fluxo " << values->FLUXO_ESCALAR[g][n] << std::endl;
            //     }
            // }

            if (aux < tolerancia)
            {
                //se a diferença entre os fluxos da iteração atual e anterior respeitar uma tolerância, então o programa para.
                break;
            }

            calculateScatteringSourceSmgi();

            iteracao++;
        }

    }
    // Cálculo final do tempo
    const auto end_time = std::chrono::high_resolution_clock::now();
    values->tempoFinalDeProcessamento = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0f;
    values->iteracaoFinal = iteracao;

    std::cout << "Iterações: " << values->iteracaoFinal << ", Tempo: " << values->tempoFinalDeProcessamento << "s"
              << ", Erro final: " << max_error << std::endl;
#else
    runDDMethodWithOneThread(data);
#endif
}

void DDMethod::runDDMethodWithTwoThreads(dados_entrada &data)
{
    /**************************************** método DD***********************************************************************

O Algoritmo

1- Colocar as condições de contorno nos fluxos angulares ψ
2- fazer a varredura para a direita e em seguida para a esquerda, obtendo fluxos angulares ψ em cada direção.
3- Calcular fluxo escalar
4- Testar convergência: a diferença entre o fluxo escalar calculado na iteração anterior com a atual deve ser menor do que
a tolerânicia para o processo parar.
5- Recalcular fonte de espalhamento e repetir todo o processo.
6- Caso a número de iterações chegue no values declarado, o processo para.

************************************************************************************************************************/
    values = &data;

    initVariables();

    const auto start_time    = std::chrono::high_resolution_clock::now();
    const long double tol    = 1.0 / pow(10, values->ordem_parada);
    int iteration            = 0;
    long double maxError     = 0.0L;

    auto rightScanFunc       = (values->tipo_ce == 2) ? &DDMethod::rightScanReflective : &DDMethod::rightScan;
    auto leftScanFunc        = (values->tipo_cd == 2) ? &DDMethod::leftScanReflective : &DDMethod::leftScan;
    auto calculateScalarFlux = (values->tipo_criterio_parada == 1) ? &DDMethod::solveScalarFluxWithRelativeCriteria
                                                                   : &DDMethod::solveScalarFluxWithAbsoluteCriteria;
    std::future<void> rightFuture;
    std::future<void> leftFuture;

    while (iteration < values->iteracao.load(std::memory_order_relaxed))
    {
        rightFuture = std::async(std::launch::async, rightScanFunc, this);
        leftFuture = std::async(std::launch::async, leftScanFunc, this);

        rightFuture.wait();
        leftFuture.wait();

        maxError = (this->*calculateScalarFlux)();

        if (maxError < tol)
        {
            break;
        }

        calculateScatteringSourceSmgi();
        iteration++;
    }

    const auto end_time = std::chrono::high_resolution_clock::now();
    values->tempoFinalDeProcessamento = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0f;
    values->iteracaoFinal = iteration;

    std::cout << "time runDDMethodWithTwoThreads() " << values->tempoFinalDeProcessamento << " criterio_parada "
              << values->tipo_criterio_parada << " maxError " << maxError
              <<" iteration "<<iteration<<std::endl;

}

void DDMethod::runDDMethodWithOneThread(dados_entrada &data)
{
/**************************************** método DD***********************************************************************

O Algoritmo

1- Colocar as condições de contorno nos fluxos angulares ψ
2- fazer a varredura para a direita e em seguida para a esquerda, obtendo fluxos angulares ψ em cada direção.
3- Calcular fluxo escalar
4- Testar convergência: a diferença entre o fluxo escalar calculado na iteração anterior com a atual deve ser menor do que
a tolerânicia para o processo parar.
5- Recalcular fonte de espalhamento e repetir todo o processo.
6- Caso a número de iterações chegue no values declarado, o processo para.

************************************************************************************************************************/
    values = &data;

    const long double tol = 1.0L / static_cast<long double>(pow(10, values->ordem_parada));
    long double maxError  = 0.0L;
    initVariables();

    auto rightScanFunc       = (values->tipo_ce == 2) ? &DDMethod::rightScanReflective : &DDMethod::rightScan;
    auto leftScanFunc        = (values->tipo_cd == 2) ? &DDMethod::leftScanReflective : &DDMethod::leftScan;
    auto calculateScalarFlux = (values->tipo_criterio_parada == 1) ? &DDMethod::solveScalarFluxWithRelativeCriteria
                                                                   : &DDMethod::solveScalarFluxWithAbsoluteCriteria;
    int iteration = 0;

    const auto start_time = std::chrono::high_resolution_clock::now();

    while (iteration < values->iteracao.load(std::memory_order_relaxed))
    {
        (this->*rightScanFunc)();
        (this->*leftScanFunc)();

        maxError = (this->*calculateScalarFlux)();

        if (maxError < tol)
        {
            break;
        }

        calculateScatteringSourceSmgi();

        iteration++;
    }

    //writeLog();

    const auto end_time               = std::chrono::high_resolution_clock::now();
    values->tempoFinalDeProcessamento = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
                                        / 1000.0f;
    values->iteracaoFinal           = iteration;

    std::cout<<"time runDDMethodWithOneThreads() "
              <<values->tempoFinalDeProcessamento
              <<" criterio_parada "<<values->tipo_criterio_parada
              <<" maxError "<<maxError
              <<" iteration "<<iteration<<std::endl;
}

void DDMethod::saveTxtFile(dados_entrada &data)
{
    values = &data;
    ofstream saida_dados;
    string titulo;
    ostringstream auxiliar;

    // Gera o nome do arquivo baseado nos parâmetros
    auxiliar << "Saida_Dados_R" << values->n_R << "_G" << values->G << "_L" << values->L << "_N" << values->n << "_Nod" << values->NODOSX << ".txt";
    titulo = auxiliar.str();

    // Abre o arquivo para escrita
    saida_dados.open(titulo);

    // Cabeçalho do arquivo
    saida_dados << string(160, '*') << endl;
    saida_dados << setw(90) << "Saída de dados" << endl;
    saida_dados << string(160, '*') << endl;

    // Informações gerais sobre o processamento
    saida_dados << "Número de iter:\t\t" << values->iteracaoFinal << endl;
    saida_dados << "Tempo(s):\t\t" << values->tempoFinalDeProcessamento / CLOCKS_PER_SEC << endl; // Converte o tempo para segundos
    saida_dados << string(160, '_') << endl;
    saida_dados << string(160, '*') << endl;

    // Título da seção de fluxo angular final
    saida_dados << "Fluxo angular final\t\titer: " << values->iteracaoFinal << endl;
    saida_dados << string(160, '*') << endl;

    // Cabeçalhos da tabela de resultados
    saida_dados << "Grupo\t\t";
    values->TAM_TOTAL = 0;

    // Calcula o comprimento total e imprime a primeira linha de distâncias
    for (int r = 0; r <= values->n_R; r++) {
        values->TAM_TOTAL += values->TAM[r];
    }

    double t = 0;
    while (t <= values->TAM_TOTAL) {
        saida_dados << setw(15)  << t << "cm";
        t += values->periodicidade;
    }
    saida_dados << endl;

    saida_dados << string(160, '*') << endl;

    // Escrita dos fluxos escalares para cada grupo
    for (int g = 0; g < values->G; g++) {
        saida_dados << g + 1 << "\t\t"; // Coluna do grupo
        int nod = 0;
        int steps = static_cast<int>(values->TAM_TOTAL / values->periodicidade);

        for (int i = 0; i <= steps; ++i)
        {
            saida_dados << setw(15) << scientific << setprecision(6) << values->FLUXO_ESCALAR[g][nod];
            nod += (values->NODOSX * values->periodicidade) / values->TAM_TOTAL;
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

    //varredura para a esquerda mi < 0
    //psi(i-1) = [(b - a) * psi(i) + Q + smgi] / (b + a)
    for (int g = 0; g < values->G; g++)
    {
        // ni representa o índice global do nodo atual na malha 1D. Neste caso, ele é decrescente
        ni = values->NODOSX;
        for (int r = (values->n_R - 1); r >= 0; r--)
        {
            Q     = values->fonte_g[g][r];
            a     = 0.5L * values->s_t[g][values->Map_R[r] - 1];
            passo = values->PASSO[r];

            for (int n = values->n_nodos[r] - 1; n >= 0; n--)
            {
                for (int o = 0; o < (values->n / 2); o++)
                {
                    b   = -values->mi[o] / passo;
                    num = ((b - a) * values->FLUXO_ANGULAR[g][ni][o]) + Q + values->smgi[g][ni - 1][o];
                    den = b + a;
                    fi  = num / den;
                    values->FLUXO_ANGULAR[g][ni - 1][o] = fi;
                }
                ni--;
            }
        }
    }
}

void DDMethod::rightScan()
{
    int ni                        = 0;
    long double fi                = 0.0L;
    long double Q                 = 0.0L;
    long double a                 = 0.0L;
    long double b                 = 0.0L;
    long double num               = 0.0L;
    long double den               = 1.0L;
    long double passo             = 0.5L;
    short int stRegionIndice      = 0;

    //varredura para a direita mi > 0
    // psi(i+1) = [(b - a) * psi(i) + Q + smgi] / (b + a)
    for (int g = 0; g < values->G; g++)
    {
        // ni representa o índice global do nodo atual na malha 1D. Neste caso, ele é crescente
        ni = 0;
        for (int r = 0; r < values->n_R; r++)
        {
            stRegionIndice = values->Map_R[r] - 1;

            Q              = values->fonte_g[g][r];
            a              = 0.5L * values->s_t[g][stRegionIndice];
            passo          = values->PASSO[r];

            for (int n = 0; n < (values->n_nodos[r]); n++)
            {
                //n muda conforme mudamos de regiao, ou seja, cada regiao tem uma quantidade de nodos
                for (int o = (values->n / 2); o < values->n; o++)
                {
                    //a varredura pela direita contempla a metadade das direções mi
                    b   = values->mi[o] / passo;
                    num = ((b - a) * values->FLUXO_ANGULAR[g][ni][o]) + Q + values->smgi[g][ni][o];
                    den = b + a;
                    fi  = num / den;
                    values->FLUXO_ANGULAR[g][ni + 1][o] = fi;
                }
                std::cout<<" noder "<<ni<<std::endl;

                ni++;
            }
        }
    }
}

void DDMethod::leftScanReflective()
{
    const double metadeQuadratura = values->n * 0.5;
    int refletido                 = 0;

    // Ex. N = 4, o = 2, a = 1 -> o = 3, a = 0
    // mi -0.86113631159405257254 w 0.34785484513745373869 0 entra
    // mi -0.33998104358485631282 w 0.65214515486254642784 1 entra
    // mi 0.33998104358485631282  w 0.65214515486254642784 2 sai
    // mi 0.86113631159405257254  w 0.34785484513745373869 3 sai
    for (int g = 0; g < values->G; g++)
    {
        refletido = values->n - 1;
        for (int o = 0; o < metadeQuadratura; o++)
        {
            //std::cout << "mi = " << values->mi[refletido] << " -m = " << values->mi[o] << std::endl;
            values->FLUXO_ANGULAR[g][values->NODOSX][o] = values->FLUXO_ANGULAR[g][values->NODOSX][refletido];
            refletido--;
        }
    }

    leftScan();
}

void DDMethod::rightScanReflective()
{
    const double metadeQuadratura = values->n * 0.5;
    int refletido                 = 0;

    // Ex. N = 4, o = 2, a = 1 -> o = 3, a = 0
    // mi -0.86113631159405257254 w 0.34785484513745373869 0 sai
    // mi -0.33998104358485631282 w 0.65214515486254642784 1 sai
    // mi 0.33998104358485631282  w 0.65214515486254642784 2 entra
    // mi 0.86113631159405257254  w 0.34785484513745373869 3 entra
    for (int g = 0; g < values->G; g++)
    {
        refletido = metadeQuadratura - 1;
        for (int o = metadeQuadratura; o < values->n; o++)
        {
            //std::cout << "mi = " << values->mi[refletido] << " -m = " << values->mi[o] << std::endl;
            values->FLUXO_ANGULAR[g][0][o] = values->FLUXO_ANGULAR[g][0][refletido];
            refletido--;
        }
    }

    rightScan();
}

void DDMethod::calculateScatteringSourceSmgi()
{
    //recalcular smgi
    long double soma2 = 0.0L;
    long double soma3 = 0.0L;

    for (int g = 0; g < values->G; g++)
    {
        // ni representa o índice global do nodo atual na malha 1D.
        int ni = 0;

        for (int r = 0; r < values->n_R; r++)
        {
            int z = values->Map_R[r] - 1;

            for (int nodo = 0; nodo < values->n_nodos[r]; nodo++)
            {
                for (int m = 0; m < values->n; m++)
                {
                    long double soma3 = 0.0L;

                    for (int g1 = 0; g1 < values->G; g1++)
                    {
                        for (int n = 0; n < values->n; n++)
                        {
                            //Fluxo interpolado deveria ser a média, mas para termos de otimização, a divisão por dois está abaixo
                            long double fluxoInterp = values->FLUXO_ANGULAR[g1][ni][n] + values->FLUXO_ANGULAR[g1][ni + 1][n];
                            long double pesoDir     = values->w[n];
                            long double termoEspalhamentoAnisotropico = 0.0L;

                            for (int l = 0; l <= values->L; l++)
                            {
                                // componenteAngularEspalhamento = soma_l( (2l+1) * P_l(mu_m) * P_l(mu_n) * sigma_s_l(g1 ->g) )
                                termoEspalhamentoAnisotropico += values->legendreProj[m][n][l] * values->s_s[g][g1][z][l];
                            }

                            soma3 += fluxoInterp * pesoDir * termoEspalhamentoAnisotropico;
                        }
                    }

                    //0.5 de Legendre vezes 0.5 da média do fluxo angular (fluxoInterp/2)
                    values->smgi[g][ni][m] = 0.25L * soma3;

                    if (std::isnan(values->smgi[g][ni][m]))
                    {
                        std::cerr << "The smgi is nan" << std::endl;
                    }
                }

                ni++;
            }
        }
    }
}

void DDMethod::writeLog()
{
    for (int g = 0; g < values->G; g++)
    {
        cout << " Grupo " << values->G << "\n\n";
        for (int n = 0; n <= values->NODOSX; n++)
        {
            cout << " nodo " << n << "fluxo " << values->FLUXO_ESCALAR[g][n] << endl;
        }
    }

    cout << " Fluxo angular " << "\n\n";

    for (int g = 0; g < values->G; g++)
    {
        cout << " Grupo " << values->G << "\n\n";

        for (int node = 0; node <= values->NODOSX; node++)
        {
            for (int n = 0; n < values->n; n++)
            {
                long double fluxo = values->FLUXO_ANGULAR[g][node][n];

                cout << " nodo " << node <<  " fluxo " << fluxo << "\n\n";
            }
        }
    }

    std::ofstream logFile;

    logFile.open("logfile.txt", std::ios::app); // Abre o arquivo em modo de adição

    if (!logFile.is_open())
    {
        std::cerr << "Erro ao abrir o arquivo de log!" << std::endl;
        return;
    }

    logFile << "Inicializando variáveis...\n";

    for (int g = 0; g < values->G; g++)
    {
        logFile << " Grupo " << g << "\n\n";
        for (int n = 0; n <= values->NODOSX; n++)
        {
            logFile << " nodo " << n << " fluxo " << values->FLUXO_ESCALAR[g][n] << std::endl;
        }
    }

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

    values->legendreProj.clear();
    values->legendreProj = std::vector<std::vector<std::vector<long double>>>(
        values->n, std::vector<std::vector<long double>>(
            values->n, std::vector<long double>(values->L + 1)));

    // Precalcula a parte (2l+1) * P_l(mu_m) * P_l(mu_n) para todas as direcoes e ordens
    for (int m = 0; m < values->n; m++)
    {
        for (int n = 0; n < values->n; n++)
        {
            for (int l = 0; l <= values->L; l++)
            {
                values->legendreProj[m][n][l] = (2 * l + 1) * values->Mat_Legendre[m][l] * values->Mat_Legendre[n][l];
            }
        }
    }
}

long double DDMethod::solveScalarFluxWithAbsoluteCriteria()
{
    long double somatorio1 = 0.0L;
    long double maxError   = -1.0L;

    for (int g = 0; g < values->G; g++)
    {
        for (int n = 0; n <= values->NODOSX; n++)
        {
            somatorio1 = 0.0L;

            for (int o = 0; o < values->n; o++)
            {
                somatorio1 = somatorio1 + (values->FLUXO_ANGULAR[g][n][o]) * (values->w[o]);
            }

            somatorio1     = somatorio1 * 0.5L;
            values->MODULO = fabs(values->FLUXO_ESCALAR[g][n] - somatorio1);
            maxError       = std::max(maxError, values->MODULO); //salva-se o maior erro entre os fluxos escalares

            if (    ( somatorio1 < 0.0L)
                || ( std::isnan( somatorio1 ) ) )
            {
                std::cerr << "Scalar flux is nan or negative" << std::endl;
                values->FLUXO_ESCALAR[g][n] = 0.0L;
            }
            else
            {
                //fluxo escalar nos nodos (não é o médio)
                values->FLUXO_ESCALAR[g][n] = somatorio1;
            }
        }
    }

    return maxError;
}

long double DDMethod::solveScalarFluxWithRelativeCriteria()
{
    long double somatorio1 = 0.0L;
    long double maxError   = -1.0L;

    for (int g = 0; g < values->G; g++)
    {
        for (int n = 0; n <= values->NODOSX; n++)
        {
            somatorio1 = 0.0L;

            for (int o = 0; o < values->n; o++)
            {
                somatorio1 = somatorio1 + (values->FLUXO_ANGULAR[g][n][o]) * (values->w[o]);
            }

            somatorio1 = somatorio1 * 0.5L;

            if (somatorio1 != 0)
            {
                values->MODULO = fabs((values->FLUXO_ESCALAR[g][n] - somatorio1) / somatorio1);
            }
            else
            {
                values->MODULO = fabs(values->FLUXO_ESCALAR[g][n] - somatorio1);
            }

            maxError = std::max(maxError, values->MODULO); //salva-se o maior erro entre os fluxos escalares

            if (    ( somatorio1 < 0.0L)
                 || ( std::isnan( somatorio1 ) ) )
            {
                std::cerr << "Scalar flux is nan or negative" << std::endl;
                values->FLUXO_ESCALAR[g][n] = 0.0L;
            }
            else
            {
                //fluxo escalar nos nodos (não é o médio)
                values->FLUXO_ESCALAR[g][n] = somatorio1;
            }
        }
    }

    return maxError;
}
