#include <iostream>
#include <fstream>         
#include <ctime>           
#include <cmath>           
#include <sstream>         
#include <vector>         
#include <iomanip>         
#include "Estruturas.h"   
#include "Legendre.h"     
#include "quadrule.hpp"    
#include "Const_matriz.h"  
#include "metodoDD.h"      

using namespace std;

int main()
{
    construir_dados("Ricardo_Barros2.txt", valor);
    DD();

    return 0;
}

