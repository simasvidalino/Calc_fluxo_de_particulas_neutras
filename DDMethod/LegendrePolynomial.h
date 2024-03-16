#pragma once

//==================================================================
/**
 *  legendre.h -- C++ functions to evaluate Legendre polynomials
 *
 *  Copyright (C) 2005 by James A. Chappell
 *
 *  Permission is hereby granted, free of charge, to any person
 *  obtaining a copy of this software and associated documentation
 *  files (the "Software"), to deal in the Software without
 *  restriction, including without limitation the rights to use,
 *  copy, modify, merge, publish, distribute, sublicense, and/or
 *  sell copies of the Software, and to permit persons to whom the
 *  Software is furnished to do so, subject to the following
 *  condition:
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *  OTHER DEALINGS IN THE SOFTWARE.
 */
//=================================================================
/*
 * legendre.h:  Version 0.01
 * Created by James A. Chappell
 * Created 29 September 2005
 *
 * History:
 * 29-sep-2005  created
 * Changed by Andréia Simas
 */
//==============

namespace Legendre
{
  inline void Pn(unsigned int n, double x, double *pn) //aqui será construído o vetor pn
  {
    if (x == 1.0){
        for (unsigned int l = 0 ; l <= n ; l++){
          pn[l] = 1.0 ;
        }
    }
    if (x == -1.0){
       for (unsigned int l = 0 ; l <= n ; l++){
          pn[l] = ((n % 2 == 0) ? 1.0 : -1.0) ;
        }
    }

    if ((x == 0.0) && (n % 2)){
      for (unsigned int l = 0 ; l <= n ; l++){
          pn[l] = 0.0 ;
        }
    }
    double pnm1=x ;
    double pnm2=1 ;
    pn[0]=1;
    pn[1]=x;

   for (unsigned int l = 2 ; l <= n ; l++){
      pn[l] = (((2.0 * (double)l) - 1.0) * x * pnm1 -
            (((double)l - 1.0) * pnm2)) / (double)l ;
      pnm2 = pnm1;
      pnm1 = pn[l];
    }
  }
}
