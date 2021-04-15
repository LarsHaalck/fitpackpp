#include "fitpackpp/BSplineCurve.h"
#include "FCMangle.h"

#include <cmath>
#include <stdexcept>

extern "C"
{
    void curfit(int* iopt, int* m, double* x, double* y, double* w, double* xb, double* xe, int* k,
        double* s, int* nest, int* n, double* t, double* c, double* fp, double* wrk, int* lwrk,
        int* iwrk, int* ier);
    void splev(
        double* t, int* n, double* c, int* k, double* x, double* y, int* m, int* e, int* ier);
}

namespace fitpackpp
{

BSplineCurve::BSplineCurve(int degree, const Vec& knots, const Vec& coeffs)
    : mDegree(degree)
    , mKnots(knots)
    , mCoeffs(coeffs)
{
}

BSplineCurve BSplineCurve::splrep(Vec x, Vec y, std::optional<Vec> w, std::optional<double> xb,
    std::optional<double> xe, int k, int task, std::optional<double> s, std::optional<Vec> t,
    bool per)
{
    int m = x.size();
    if (!w)
    {
        w = Vec(m, 1.0);
        if (!s)
            s = 0.0;
    }
    else if (!s)
        s = m - std::sqrt(2 * m);

    // check if input is valid
    if (static_cast<int>(w.value().size()) != m || static_cast<int>(y.size()) != m)
        throw std::runtime_error("sizes of weight, x and y musts be equal");
    if (m <= k)
        throw std::runtime_error("m > k must hold");
    if (task > 0)
        throw std::runtime_error("Task must be -1, 0");

    if (!xb)
        xb = x.front();
    if (!xe)
        xe = x.back();

    int nest = 0;
    if (task == -1)
    {
        if (!t)
            throw std::runtime_error("Knots must be given for task=-1");

        auto numKnots = t.value().size();
        Vec newT(numKnots + 2 * k + 2);
        std::copy(std::begin(t.value()), std::end(t.value()), std::begin(newT) + k + 1);
        t.value().swap(newT);
        nest = t.value().size();
    }
    else if (task == 0)
    {
        if (per)
            nest = std::max(m + 2 * k, 2 * k + 3);
        else
            nest = std::max(m + k + 1, 2 * k + 3);
        t = Vec(nest);
    }

    Vec wrk;
    if (per)
        wrk = Vec(m * (k + 1) + nest * (8 + 5 * k));
    else
        wrk = Vec(m * (k + 1) + nest * (7 + 3 * k));

    auto iwrk = std::vector<int>(nest);
    int lwrk = wrk.size();

    auto coeffs = Vec(nest);
    double fp = 0.0;
    int n = 0;

    int ier = 0;
    curfit(&task, &m, x.data(), y.data(), w.value().data(), &xb.value(), &xe.value(), &k,
        &s.value(), &nest, &n, t.value().data(), coeffs.data(), &fp, wrk.data(), &lwrk, iwrk.data(),
        &ier);

    // truncate if needed
    t.value().resize(n);
    coeffs.resize(n);
    return BSplineCurve(k, t.value(), coeffs);
}

Vec BSplineCurve::evaluate(Vec x, int ext)
{
    Vec y(x.size(), 0.0);
    int n = mKnots.size();
    int m = x.size();
    int ier = 0;

    splev(mKnots.data(), &n, mCoeffs.data(), &mDegree, x.data(), y.data(), &m, &ext, &ier);

    if (ier > 0)
    {
        if (ier == 1)
            throw std::runtime_error("Argument out of bounds");
        if (ier == 10)
            throw std::runtime_error("Invalid input data");
    }
    return y;
}
} // namespace fitpackpp
