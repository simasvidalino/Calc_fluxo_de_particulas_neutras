#pragma once

namespace Legendre
{
inline void Pn(unsigned int n, double x, double* pn)
{
    if (n == 0) {
        pn[0] = 1.0;
        return;
    }

    if (x == 1.0) {
        for (unsigned int l = 0; l <= n; ++l)
            pn[l] = 1.0;
        return;
    }

    if (x == -1.0) {
        for (unsigned int l = 0; l <= n; ++l)
            pn[l] = (l % 2 == 0 ? 1.0 : -1.0);
        return;
    }

    pn[0] = 1.0;
    if (n >= 1) pn[1] = x;

    for (unsigned int l = 2; l <= n; ++l)
        pn[l] = (((2.0 * l - 1.0) * x * pn[l - 1]) -
                 ((l - 1.0) * pn[l - 2])) / l;
}
}
