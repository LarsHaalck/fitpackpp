#ifndef BSPLINE_CURVE_H
#define BSPLINE_CURVE_H

#include <optional>
#include <vector>

namespace fitpackpp
{
using Vec = std::vector<double>;
class BSplineCurve
{
public:
    BSplineCurve(int degree, const Vec& knotX, const Vec& coefs);

    static BSplineCurve splrep(Vec x, Vec y, std::optional<Vec> w = {},
        std::optional<double> xb = {}, std::optional<double> xe = {}, int k = 3, int task = 0,
        std::optional<double> s = {}, std::optional<Vec> t = {}, bool per = false);

    Vec evaluate(Vec x, int e = 0);
    Vec operator()(Vec x, int e = 0) { return evaluate(x, e); };
    double evaluate(double x, int e = 0) { return evaluate(Vec {x}, e)[0]; };
    double operator()(double x, int e = 0) { return evaluate(Vec {x}, e)[0]; };

    int degree() const { return mDegree; };
    std::vector<double> knots() const { return mKnots; }
    std::vector<double> coeffs() const { return mCoeffs; }

private:
    int mDegree;
    std::vector<double> mKnots;
    std::vector<double> mCoeffs;
};

} // namespace fitpackpp

#endif // BSPLINE_CURVE_H
