#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>

#define PI 3.14159265358979323846

std::tuple<std::vector<double>, std::vector<double>> get_GQ(int N)
{
    std::vector<double> mi(N);
    std::vector<double> w(N);

    int M = (N + 1) / 2;
    double u, u1, P1, P2, P3, DP;

    for (int i = 0; i < M; i++) {
        u = std::cos(PI * (i + 0.75) / (N + 0.5));

        do {
            P1 = 1.0;
            P2 = 0.0;
            for (int j = 0; j < N; j++) {
                P3 = P2;
                P2 = P1;
                P1 = ((2.0 * j + 1.0) * u * P2 - j * P3) / (j + 1.0);
            }

            DP = N * (u * P1 - P2) / (u * u - 1.0);
            u1 = u;
            u = u1 - P1 / DP;
        } while (std::fabs(u - u1) > 1e-15);

        double weight = 2.0 / ((1.0 - u * u) * DP * DP);

        mi[i]         = -u;
        mi[N - i - 1] = u;
        w[i]          = weight;
        w[N - i - 1]  = weight;
    }

    return {mi, w};
}
