#ifndef __FUNCTOR_CLASS_H__
#define __FUNCTOR_CLASS_H__

#include <eigen3/Eigen/Eigenvalues>
#include "forceFunction.hpp"
#include <boost/numeric/odeint.hpp>
#include "Helpers.hpp"

class functorClass
{
    private:
        /* @brief Returns the sign of a given input variable
           @param[in] val numeric type to return the sign of
           @returns {+1, 0, -1}: Positive, zero, negative, respectively.
        */
        template <typename numericType>
        int sgn(numericType val)
        {
            return (val > 0) - (val < 0);
        }

        /* Integrate a given std::vector<double> forward from initialTime to finalTime.
        @param[inout] x std::vector<double> containing the state to integrate. On entry, it is the initial condition. On exit, it is the final condition.
        */
        void integrate (std::vector<double>& x)
        {
            typedef std::vector<double> state_type;
            int sign = sgn(this->finalTime - this->initialTime);
            boost::numeric::odeint::bulirsch_stoer<state_type> bulirsch_stepper(this->absTol, this->relTol);
            boost::numeric::odeint::integrate_adaptive(bulirsch_stepper, dynSystem, x, this->initialTime, this->finalTime, sign*.01, abcFlowObserver);
        }

    public:
        Eigen::Vector3d previousSolution;   // For globally smooth strainlines
        double initialTime;                 // For numerical integration
        double finalTime;                   // ""
        double xGridSpacing;                // For computing derivatives
        double yGridSpacing;                // ""
        double zGridSpacing;                // ""
        double absTol = 1e-012;             // Absolute integration tolerance
        double relTol = 1e-012;             // Relative integration tolerance

    /* Class constructor */
    functorClass(){}; // Blank

    /* @brief Force function for repelling strainline advection.
       @param[in] x The value of the state at a given time.
       @param[inout] xdot The value of the state derivative at a given time.
       @param[in] t The current time of the integration.
    */
    void operator()(const std::vector<double>& x, std::vector<double>& xdot, double t)
    {
        /* We need to know the direction of the eigenvector field at this point.
         *
         * Therefore, construct a grid around this initial point and integrate the trajectories forward
         */
        std::vector<double> left(3), up(3), right(3), down(3), pZ(3), mZ(3);
        Eigen::Vector3d direction;

        double maxEigenvalue, minEigenvalue;

        /* Construct initial conditions grid */
        left = x; left[0] -= xGridSpacing;
        right = x; right[0] += xGridSpacing;
        up = x; up[1] += yGridSpacing;
        down = x; down[1] -= yGridSpacing;
        pZ = x; pZ[2] += zGridSpacing;
        mZ = x; mZ[2] -= zGridSpacing;

        /* Construct the normal to the plane */
        Eigen::Vector3d normal, unitNormal, xplane, yplane;
        for (size_t i = 0; i < left.size(); ++i)
        {
            xplane[i] = right[i] - left[i];
            yplane[i] = up[i] - down[i];
        }
        normal = xplane.cross(yplane);
        unitNormal = normal.normalized();

        /* Integrate the initial conditions forward */
        this->integrate(left);
        this->integrate(right);
        this->integrate(up);
        this->integrate(down);
        this->integrate(pZ);
        this->integrate(mZ);

        /* Compute the dominant eigenvector */
        Helpers::getDominantEigenvectorAndEigenvalue(left, up, right, down, pZ, mZ, xGridSpacing, yGridSpacing, zGridSpacing, 
                                                     direction, maxEigenvalue, minEigenvalue);

        /* Construct the derivative vector.
         * Recall that the equations of motion for the strainline are \gamma^\prime = n_\Pi \cross \eta_3 
         */
        Eigen::Vector3d derivative = unitNormal.cross(direction);
        double innerProduct = derivative.dot(this->previousSolution);
        derivative = derivative * sgn(innerProduct);
        /* Assign */
        xdot[0] = derivative[0]; 
        xdot[1] = derivative[1];
        xdot[2] = derivative[2];
    }
};

#endif