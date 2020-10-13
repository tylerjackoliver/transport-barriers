#ifndef __HELPER_ROUTINES_H__
#define __HELPER_ROUTINES_H__

#include<eigen3/Eigen/Eigenvalues>
#include<boost/numeric/odeint.hpp>

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

}

#endif