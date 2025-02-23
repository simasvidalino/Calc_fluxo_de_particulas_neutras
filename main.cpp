#include <string>
#include "BuildMatricesWIthTxt.h"//struct para guardar variáveis
#include "VariablesUsed.h"
#include "DDNumericalMethod.h"

using namespace std;

int main()
{
    static dados_entrada valor;
    string path = "/home/andreia-simas/Documentos/TCC//Arquivo_de_entrada.txt";
    char continuar;

    std::cout<<"Informe o caminho completo do arquivo de entrada"<<std::endl;
    std::cin>>path;

    do {
        std::cout << "Calcular o Fluxo médio de partículas neutras: ";

        try
        {
            construir_dados(path, valor);
            DDMethod::getInstance()->runDDMethodWithOneThread(valor);
            DDMethod::getInstance()->saveTxtFile(valor);

            std::cout << "Deseja calcular novamente? (s/n): ";
            std::cin >> continuar;
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << "Exception caught: " << e.what() << '\n';
        }

    } while (continuar == 's' || continuar == 'S');



    return 0;
}
