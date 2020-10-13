#ifndef __FUNCTOR_CLASS_H__
#define __FUNCTOR_CLASS_H__

#include <eigen3/Eigen/Eigenvalues>
#include "forceFunction.hpp"
#include <boost/numeric/odeint.hpp>
#include "Helpers.hpp"


class functorClass
{
    private:
        /* Return the sign of a given input x; return [1, 0, -1] */
        template <typename timeType>
        int sgn(timeType val)
        {
            return (val > 0) - (val < 0);
        }

        /* Integrate a given state x forward from initialTime to finalTime */
        // template <typename Type>
        void integrate (std::vector<double>& x)
        {
            typedef std::vector<double> state_type;
            int sign = sgn(this->finalTime - this->initialTime);
            // boost::numeric::odeint::bulirsch_stoer<state_type> bulirsch_stepper(this->absTol, this->relTol);
            // boost::numeric::odeint::integrate_adaptive(bulirsch_stepper, dynSystem, x, this->initialTime, this->finalTime, sign*.01, abcFlowObserver);
        }

    public:
        Eigen::Vector3d previousSolution;   // For globally smooth strainlines
        double initialTime;               // For numerical integration
        double finalTime;                 // ""
        double xGridSpacing;                // For computing derivatives
        double yGridSpacing;                // ""
        double zGridSpacing;                // ""
        double absTol = 1e-012;             // Absolute integration tolerance
        double relTol = 1e-012;             // Relative integration tolerance

    /* Class constructor */
    functorClass(){}; // Blank

    /* requires arguments to be vectors */
    template <typename Type>
    void operator()(std::vector<Type>& x, std::vector<Type>& xdot, Type t)
    {
        /* We need to know the direction of the eigenvector field at this point - therefore, construct a grid around this initial point and integrate the trajectories forward */
        std::vector<double> left(3), up(3), right(3), down(3), pZ(3), mZ(3);
        Eigen::Vector3d direction;

        left = x; left[0] -= xGridSpacing;
        right = x; right[0] += xGridSpacing;
        up = x; up[1] += yGridSpacing;
        down = x; down[1] -= yGridSpacing;
        pZ = x; pZ[2] += zGridSpacing;
        mZ = x; mZ[2] -= zGridSpacing;

        // Integrate these all forwards
        this->integrate(left);
        this->integrate(right);
        this->integrate(up);
        this->integrate(down);
        this->integrate(pZ);
        this->integrate(mZ);

        // Compute the eigenvector
        Helpers::getDominantEigenVector(left, up, right, down, pZ, mZ, xGridSpacing, yGridSpacing, zGridSpacing, direction);

        // Create the derivative vector
        Type innerProduct = this->previousSolution.dot(direction);
        Eigen::Vector3d derivative = innerProduct * direction;
        xdot[0] = derivative[0]; 
        xdot[1] = derivative[1];
        xdot[2] = derivative[2];
    }
};

#endif