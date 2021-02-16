#ifndef __FORCE_FUNCTION_H__
#define __FORCE_FUNCTION_H__

#include <vector>
#include <cmath>

typedef std::vector<double> state_type;

void abcFlow(state_type& x, std::vector<double>& xdot, const double t)
{
    double A = std::sqrt(3);
    double B = std::sqrt(2);
    double C = 1.0;

    xdot[0] = A * std::sin(x[2]) + C * std::cos(x[1]);
    xdot[1] = B * std::sin(x[0]) + A * std::cos(x[2]);
    xdot[2] = C * std::sin(x[1]) + B * std::cos(x[0]);
};

void pdg(state_type& x, std::vector<double>& xdot, const double t)
{
    double epsilon = 0.25;
    double A = 0.1;
    double omega = 4.0 * std::atan(1.0) / 10.0;

    double at = epsilon * std::sin(omega * t);
    double bt = 1 - 2 * epsilon * std::sin(omega * t);
    double f = at*x[0]*x[0] + bt*x[0];
    double dfdx = 2*at*x[0] + bt;

    double pi = 4 * std::atan(1);
    double u = -pi * A * std::sin(pi*f) * std::cos(pi*x[1]);
    double v = pi * A * std::cos(pi*f) * std::sin(pi*x[1]) * dfdx;

    xdot[0] = u;
    xdot[1] = v;
    xdot[2] = 0.0;
};

void periodicABCFlow(state_type& x, state_type& xdot, const double t)
{
    double A = std::sqrt(3);
    double B = std::sqrt(2);
    double C = 1.0;

    xdot[0] = (A + 0.1 * std::sin(t) ) * std::sin(x[2]) + C * std::cos(x[1]);
    xdot[1] = B * std::sin(x[0]) + (A + 0.1 * std::sin(t) ) * std::cos(x[2]);
    xdot[2] = C * std::sin(x[1]) + B * std::cos(x[0]);
}

void dynSystem(state_type& x, std::vector<double>& xdot, const double t)
{
    periodicABCFlow(x, xdot, t);
}

// template <typename Type>
void abcFlowObserver(state_type& x, const double t){};

void doubleGyre(state_type& x, state_type& dxdt, const double t){

    double epsilon = 0.;
    double A = 0.1;
    double omega = 4.0 * std::atan(1.0) / 10.0;

    double at = epsilon * std::sin(omega * t);
    double bt = 1 - 2 * epsilon * std::sin(omega * t);
    double f = at*x[0]*x[0] + bt*x[0];
    double dfdx = 2*at*x[0] + bt;

    double pi = 4 * std::atan(1);
    double u = -pi * A * std::sin(pi*f) * std::cos(pi*x[1]);
    double v = pi * A * std::cos(pi*f) * std::sin(pi*x[1]) * dfdx;

    dxdt[0] = u;
    dxdt[1] = v;

}

#endif
