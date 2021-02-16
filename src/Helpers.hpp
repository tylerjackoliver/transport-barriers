#ifndef __HELPER_ROUTINES_H__
#define __HELPER_ROUTINES_H__

#include<eigen3/Eigen/Eigenvalues>
#include<boost/numeric/odeint.hpp>
#include<vector>
#include <numeric>
#include "forceFunction.hpp"

namespace Helpers
{
    /* @brief Computes the sign of a given number.
       @param[in] val The number to compute the sign of
       @returns Integer in {+1, 0, -1} giving the sign of val.
    */
    template <typename numericType>
    int sgn(numericType val)
    {
        return (val > 0) - (val < 0);
    }

    /* @brief Integrates a given trajectory from time initTime to finalTime
       @param[inout] x On entry, it is the initial condition. On exit, it is the final condition.
       @param[in] initTime The initial time of the integration
       @param[in] finalTime The final time of the integration.
       @param[in] absTol The absolute tolerance for the integrator
       @param[in] relTol The relative tolerance for the integrator
    */
    template <typename numericType>
    void integrateForceFunction(std::vector<numericType>& x, numericType& initTime, numericType& finalTime, double& absTol, double& relTol)
    {
        typedef std::vector<numericType> state_type;
        int sign = Helpers::sgn(finalTime - initTime);
        boost::numeric::odeint::bulirsch_stoer<state_type> bulirsch_stepper(absTol, relTol);
        boost::numeric::odeint::integrate_adaptive(bulirsch_stepper, dynSystem, x, initTime, finalTime, sign*.01, abcFlowObserver);
    }


    template <typename vecType, typename scalarType>
    void getDominantEigenVector(vecType& left, vecType& up, vecType& right, vecType& down, vecType &pZ, vecType &mZ, scalarType &xGridSpacing, scalarType &yGridSpacing, scalarType &zGridSpacing, Eigen::Vector3d& ev)
    {
        /* Construct CGST */
        Eigen::Matrix3d deformation, cauchy_green;

        // Deformation tensor
        deformation(0,0) = (right[0]-left[0]) / (2. * xGridSpacing);
        deformation(0,1) = (up[0]-down[0]) / (2. * yGridSpacing);
        deformation(0,2) = (pZ[0]-mZ[0]) / (2. * zGridSpacing);
        
        deformation(1,0) = (right[1]-left[1]) / (2. * xGridSpacing);
        deformation(1,1) = (up[1]-down[1]) / (2. * yGridSpacing);
        deformation(1,2) = (pZ[1]-mZ[1]) / (2. * zGridSpacing);

        deformation(2,0) = (right[2]-left[2]) / (2. * xGridSpacing);
        deformation(2,1) = (up[2]-down[2]) / (2. * yGridSpacing);
        deformation(2,2) = (pZ[2]-mZ[2]) / (2. * zGridSpacing);
                            
        cauchy_green = deformation.transpose() * deformation;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(cauchy_green);
        ev = solver.eigenvectors().col(2).real();
    }

    template <typename vecType, typename scalarType>
    void getDominantEigenvectorAndEigenvalue(vecType& left, vecType& up, vecType& right, vecType& down, vecType &pZ, vecType &mZ, scalarType &xGridSpacing, scalarType &yGridSpacing, scalarType &zGridSpacing, Eigen::Vector3d& ev, double& dominantEigenvalue, double& minimumEigenvalue)
    {
        /* Construct CGST */
        Eigen::Matrix3d deformation, cauchy_green;

        // Deformation tensor
        deformation(0,0) = (right[0]-left[0]) / (2. * xGridSpacing);
        deformation(0,1) = (up[0]-down[0]) / (2. * yGridSpacing);
        deformation(0,2) = (pZ[0]-mZ[0]) / (2. * zGridSpacing);
        
        deformation(1,0) = (right[1]-left[1]) / (2. * xGridSpacing);
        deformation(1,1) = (up[1]-down[1]) / (2. * yGridSpacing);
        deformation(1,2) = (pZ[1]-mZ[1]) / (2. * zGridSpacing);

        deformation(2,0) = (right[2]-left[2]) / (2. * xGridSpacing);
        deformation(2,1) = (up[2]-down[2]) / (2. * yGridSpacing);
        deformation(2,2) = (pZ[2]-mZ[2]) / (2. * zGridSpacing);

        cauchy_green = deformation.transpose() * deformation;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(cauchy_green);
        ev = solver.eigenvectors().col(2).real();
        dominantEigenvalue = solver.eigenvalues().maxCoeff();
        minimumEigenvalue = solver.eigenvalues().minCoeff();
    }

    /* @brief Get the eigenvector about a central point based on auxiliary grid
       @param[in] Point<vecType> The current point about which to get the eigenvector about
       @param[in] xStep Numeric type representing the step in the x-direction to take
       @param[in] yStep Numeric type representing the step in the y-direction to take
       @param[in] zStep Numeric type representing the step in the z-direction to take
       @param[out] vector Eigen::Matrix<double, 3, 1> representing the output Eigenvector
       @param[in] absTol The absolute tolerance for the numerical integration
       @param[in] relTol The relative tolerance for the numerical integration
       @parma[in] initTime The initial time for the integration
       @param[in] finalTime The final time for the integration
    */
    template <typename vecType, typename scalarType>
    void makeGridGetEigenVector(Point<vecType>& thisPoint, scalarType& xStep, scalarType& yStep, scalarType& zStep, Eigen::Matrix<double, 3, 1>& vector, double absTol, double relTol, double initTime, double finalTime)
    {
        /* Form the grid about thisPoint */
        Point<vecType> xLeft, xRight, yPlus, yMinus, zPlus, zMinus;
        xLeft = thisPoint; xLeft.x -= xStep;
        xRight = thisPoint; xRight.x += xStep;
        yPlus = thisPoint; yPlus.y += yStep;
        yMinus = thisPoint; yMinus.y -= yStep;
        zPlus = thisPoint; zPlus.z += zStep;
        zMinus = thisPoint; zMinus.z -= zStep;

        /* Convert to vectors */
        std::vector<double> xL(3), xR(3), yP(3), yM(3), zP(3), zM(3);
        xL[0] = xLeft.x; xL[1] = xLeft.y; xL[2] = xLeft.z;
        xR[0] = xRight.x; xR[1] = xRight.y; xR[2] = xRight.z;
        yP[0] = yPlus.x; yP[1] = yPlus.y; yP[2] = yPlus.z;
        yM[0] = yMinus.x; yM[1] = yMinus.y; yM[2] = yMinus.z;
        zP[0] = zPlus.x; zP[1] = zPlus.y; zP[2] = zPlus.z;
        zM[0] = zMinus.x; zM[1] = zMinus.y; zM[2] = zMinus.z;

        /* Integrate forward */
        integrateForceFunction(xL, initTime, finalTime, absTol, relTol);
        integrateForceFunction(xR, initTime, finalTime, absTol, relTol);
        integrateForceFunction(yP, initTime, finalTime, absTol, relTol);
        integrateForceFunction(yM, initTime, finalTime, absTol, relTol);
        integrateForceFunction(zP, initTime, finalTime, absTol, relTol);
        integrateForceFunction(zM, initTime, finalTime, absTol, relTol);

        /* From this, compute the eigenvector */
        getDominantEigenVector(xL, yP, xR, yM, zP, zM, xStep, yStep, zStep, vector);
    }

    /* @brief Gets difference between a specified point and a vector
     * @param[in] Vector<Vector>> of the strainline trajectory
     * @param[in] Vector of the desired point
     * @param[in] Distance threshold
     * @param[out] boolean whether the point is below the given distance
     */
    template <typename Type>
    bool distanceBelowThreshold(std::vector<Point<Type>>& trajectory, Point<Type>& point, Type distanceThreshold)
    {
        double distanceThresholdSquared = distanceThreshold * distanceThreshold;
        for (size_t i = 0; i < trajectory.size(); ++i)
        {
            Point<Type> differencePoint = trajectory[i] - point;
            Type distance = (differencePoint.x * differencePoint.x) + (differencePoint.y * differencePoint.y) + (differencePoint.z * differencePoint.z);
            if (distance < distanceThresholdSquared)
            {
                return true;
            }
        }
        return false;
    }

}

#endif