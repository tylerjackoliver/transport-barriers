#ifndef __HELPER_ROUTINES_H__
#define __HELPER_ROUTINES_H__

#include<eigen3/Eigen/Eigenvalues>
#include<boost/numeric/odeint.hpp>
#include<vector>
#include "forceFunction.hpp"

namespace Helpers
{

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
    
    /* @brief Get the dominant eigenvector for a given position by constructing the grid around it.
     * @param[in] x0 std::vector<T> containing the initial state to determine the eigenvector for
     * @param[inout] ev The dominant eigenvector for the flow at that point.
     * @param[in] xSpace The step in the x-direction
     * @param[in] ySpace The step in the y-direction
     * @param[in] zSpace the step in the z-direction.
     */
    template <typename Type>
    void getDominantEigenVector(const std::vector<Type>& x0, Type& xSpace, Type& ySpace, Type& zSpace, Eigen::Matrix<Type, 3, 1>& ev)
    {
        Eigen::Matrix<Type, 3, 3> deformation, cauchy_green;
        std::vector<double> xPlus(3), xMinus(3), yPlus(3), yMinus(3), zPlus(3), zMinus(3);

        xPlus = x0; xPlus[0] += xSpace;
        xMinus = x0; xMinus[0] -= xSpace;
        yPlus = x0; yPlus[1] += ySpace;
        yMinus = x0; yMinus[1] -= ySpace;
        zPlus = x0; zPlus[2] += zSpace;
        zMinus = x0; zMinus[2] -= zSpace;

        this->integrate(xPlus);
        this->integrate(xMinus);
        this->integrate(yPlus);
        this->integrate(yMinus);
        this->integrate(zPlus);
        this->integrate(zMinus);

        // deformation tensor 
        deformation(0,0) = (xPlus[0]-xMinus[0]) / (2. * xSpace);
        deformation(0,1) = (yPlus[0]-yMinus[0]) / (2. * ySpace);
        deformation(0,2) = (zPlus[0]-zMinus[0]) / (2. * zSpace);
        
        deformation(1,0) = (xPlus[1]-xMinus[1]) / (2. * xSpace);
        deformation(1,1) = (yPlus[1]-yMinus[1]) / (2. * ySpace);
        deformation(1,2) = (zPlus[1]-zMinus[1]) / (2. * zSpace);

        deformation(2,0) = (xPlus[2]-xMinus[2]) / (2. * xSpace);
        deformation(2,1) = (yPlus[2]-yMinus[2]) / (2. * ySpace);
        deformation(2,2) = (zPlus[2]-zMinus[2]) / (2. * zSpace);
        
        cauchy_green = deformation.transpose() * deformation;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(cauchy_green);
        ev = solver.eigenvectors().col(2).real();
    }

    /* @brief Get the dominant eigenvector and eigenvalue for a given position by constructing the grid around it.
     * @param[in] x0 std::vector<T> containing the initial state to determine the eigenvector for
     * @param[inout] ev The dominant eigenvector for the flow at that point.
     * @param[inout] eVal The dominant eigenvalue for the flow at that point.
     * @param[in] xSpace The step in the x-direction
     * @param[in] ySpace The step in the y-direction
     * @param[in] zSpace the step in the z-direction.
     */
    template <typename Type>
    void getDominantEigenVectorAndEigenvalue(const std::vector<Type>& x0, Type& xSpace, Type& ySpace, Type& zSpace, Eigen::Matrix<Type, 3, 1>& ev, Type& eVal)
    {
        Eigen::Matrix<Type, 3, 3> deformation, cauchy_green;
        std::vector<double> xPlus(3), xMinus(3), yPlus(3), yMinus(3), zPlus(3), zMinus(3);

        xPlus = x0; xPlus[0] += xSpace;
        xMinus = x0; xMinus[0] -= xSpace;
        yPlus = x0; yPlus[1] += ySpace;
        yMinus = x0; yMinus[1] -= ySpace;
        zPlus = x0; zPlus[2] += zSpace;
        zMinus = x0; zMinus[2] -= zSpace;

        this->integrate(xPlus);
        this->integrate(xMinus);
        this->integrate(yPlus);
        this->integrate(yMinus);
        this->integrate(zPlus);
        this->integrate(zMinus);

        // deformation tensor 
        deformation(0,0) = (xPlus[0]-xMinus[0]) / (2. * xSpace);
        deformation(0,1) = (yPlus[0]-yMinus[0]) / (2. * ySpace);
        deformation(0,2) = (zPlus[0]-zMinus[0]) / (2. * zSpace);
        
        deformation(1,0) = (xPlus[1]-xMinus[1]) / (2. * xSpace);
        deformation(1,1) = (yPlus[1]-yMinus[1]) / (2. * ySpace);
        deformation(1,2) = (zPlus[1]-zMinus[1]) / (2. * zSpace);

        deformation(2,0) = (xPlus[2]-xMinus[2]) / (2. * xSpace);
        deformation(2,1) = (yPlus[2]-yMinus[2]) / (2. * ySpace);
        deformation(2,2) = (zPlus[2]-zMinus[2]) / (2. * zSpace);
        
        cauchy_green = deformation.transpose() * deformation;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(cauchy_green);
        ev = solver.eigenvectors().col(2).real();
        eVal = solver.eigenvalues().maxCoeff();
    }

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

    /* Construct helicity by construction obtaining directly the eigenvectors around the desired points using an auxiliary grid. */
    template <typename Type>
    void getHelicityAuxiliaryGrid(std::vector<Type>& x, Type& initTime, Type& finalTime, double& absTol, double&relTol,
                                  double &xStep, double& yStep, double& zStep)
    {
        /* Define auxiliary grid of points */

    }

}

#endif