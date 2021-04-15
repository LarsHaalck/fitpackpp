#include <algorithm>
#include <iostream>
#include <fstream>

#include "fitpackpp/BSplineCurve.h"

using namespace fitpackpp;

Vec linspace(double a, double b, std::size_t n)
{
    double delta = (b - a) / static_cast<double>(n - 1);
    std::vector<double> values(n);
    double val = a;
    for (auto it = std::begin(values); it != std::end(values); ++it)
    {
        *it = val;
        val += delta;
    }
    return values;
}

int main()
{
    Vec x {1, 2, 3, 4, 5};
    Vec y {4, 0, 1, 3, 6};
    Vec w {1, 1, 0, 1, 2};

    auto spline = BSplineCurve::splrep(x, y, w, x.front(), x.back());

    for (auto e : spline.knots())
        std::cout << e << ", ";
    std::cout << std::endl;

    for (auto e : spline.coeffs())
        std::cout << e << ", ";
    std::cout << std::endl;

    std::cout << spline.degree() << std::endl;

    auto xx = linspace(0, 6, 200);
    auto yy = spline(xx);

    std::ofstream csv("spline.csv");
    for (std::size_t i = 0; i < yy.size(); i++)
        csv << xx[i] << "," << yy[i] << "\n";

}
