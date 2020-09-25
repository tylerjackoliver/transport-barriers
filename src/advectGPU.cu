/*
 Copyright 2011-2012 Karsten Ahnert
 Copyright 2011-2013 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <iostream>
#include <cmath>
#include <utility>
#include <chrono>

#include <cuda.h>
#include <omp.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

#include <boost/numeric/odeint.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>

struct xyz
{

    double x;
    double y;
    double z;

};

// CUDA definitions
typedef double value_type;
typedef thrust::host_vector<value_type> state_typeCPU;
typedef thrust::device_vector<value_type> state_typeGPU;
typedef boost::multi_array<xyz, 3> field_;

struct abcSystem{

    struct abcFunctor
    {

        template<class T>
        __host__ __device__ 
        void operator()(T t) const{

            // Unpack the parameters
            value_type x = thrust::get<0>(t);
            value_type y = thrust::get<1>(t);
            value_type z = thrust::get<2>(t);

            thrust::get<3>(t) = sqrt(3.0) * sin(z) + 1.0 * cos(y);
            thrust::get<4>(t) = sqrt(2.0) * sin(x) + sqrt(3.0) * cos(z);
            thrust::get<5>(t) = 1.0 * sin(y) + sqrt(2) * cos(x);

        }

    };

    abcSystem(size_t N) : dim_N(N){}

    template <class state, class deriv>
    void operator()(const state& x, deriv& xdot, value_type t) const
    {

        thrust::for_each(

            thrust::make_zip_iterator( thrust::make_tuple(

                boost::begin(x),
                boost::begin(x) + dim_N,
                boost::begin(x) + 2 * dim_N,
                boost::begin(xdot),
                boost::begin(xdot) + dim_N,
                boost::begin(xdot) + 2 * dim_N)
                
            ),

            thrust::make_zip_iterator(thrust::make_tuple(

                boost::begin(x) + dim_N,
                boost::begin(x) + 2 * dim_N,
                boost::begin(x) + 3 * dim_N,
                boost::begin(xdot) + dim_N,
                boost::begin(xdot) + 2 * dim_N,
                boost::begin(xdot) + 3 * dim_N )

            ), abcFunctor() );

    };

    size_t dim_N;

};
 
const value_type sigma = 10.0;
const value_type b = 8.0 / 3.0;

//[ thrust_lorenz_parameters_define_simple_system
struct lorenz_system
{
    struct lorenz_functor
    {
        template< class T >
        __host__ __device__
        void operator()( T t ) const
        {
            // unpack the parameter we want to vary and the Lorenz variables
            value_type x = thrust::get< 0 >( t );
            value_type y = thrust::get< 1 >( t );
            value_type z = thrust::get< 2 >( t );
            thrust::get< 4 >( t ) = 1.0 * ( y - x );
            thrust::get< 5 >( t ) = x - y - x * z;
            thrust::get< 6 >( t ) = -1.0 * z + x * y ;

        }
    };

    lorenz_system( size_t N )
    : m_N( N ){ }

    template< class State , class Deriv >
    void operator()(  const State &x , Deriv &dxdt , value_type t ) const
    {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple(
                        boost::begin( x ) ,
                        boost::begin( x ) + m_N ,
                        boost::begin( x ) + 2 * m_N ,
                        boost::begin( dxdt ) ,
                        boost::begin( dxdt ) + m_N ,
                        boost::begin( dxdt ) + 2 * m_N  ) ) ,
                thrust::make_zip_iterator( thrust::make_tuple(
                        boost::begin( x ) + m_N ,
                        boost::begin( x ) + 2 * m_N ,
                        boost::begin( x ) + 3 * m_N ,
                        boost::begin( dxdt ) + m_N ,
                        boost::begin( dxdt ) + 2 * m_N ,
                        boost::begin( dxdt ) + 3 * m_N  ) ) ,
                lorenz_functor() );
    }
    size_t m_N;
};
 
void advectPositionGPUDriver(field_& initPos, field_& finalPos, double initTime, double finalTime, \
double absTol, double relTol)
 {
    int driverVersion, runtimeVersion;
    cudaDriverGetVersion(&driverVersion); cudaRuntimeGetVersion(&runtimeVersion);

    std::cout << "Computing on CUDA version \t" << driverVersion << "\t and runtime version \t" << runtimeVersion << std::endl;

    // Get problem dimensionality
    unsigned int xSize = initPos.shape()[0];
    unsigned int ySize = initPos.shape()[1];
    unsigned int zSize = initPos.shape()[2];
    // Get the size of the flattened array
    unsigned int extent = xSize * ySize * zSize;
    // Number of dimensions for the array
    unsigned int dim = 3;
    // Assign flattened array
    state_typeCPU x(dim*extent);

    std::cout << "The total number of values is " << extent << std::endl;

    /* Fill the device vector with the correct values.
    *  For speed, initialise _first_ on the host device (x), then transfer to the GPU in one chunk (dx)
    */

    std::cout << "Transferring data into CPU array...";
    auto start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for (int i = 0; i < (int)xSize-1; i++)
    {  
        for (int j = 0; j < (int)ySize-1; j++)
        {
            for (int k = 0; k < (int)zSize-1; k++)
            {

                xyz tmp = initPos[i][j][k];

                x[i * xSize * ySize + j * ySize + k] = tmp.x;
                x[extent + i * xSize * ySize + j * ySize + k] = tmp.y;
                x[2*extent + i * xSize * ySize + j * ySize + k] = tmp.z; 

            } // k

        } // j

    } // i

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "done. Time required: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms." << std::endl;

    // Now copy data to GPU array

    std::cout << "Copying data from host to device...";
    start = std::chrono::high_resolution_clock::now();
    state_typeGPU dX = x;
    end = std::chrono::high_resolution_clock::now();
    std::cout << "done. Time required: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms." << std::endl;

    boost::numeric::odeint::bulirsch_stoer<state_typeGPU> stepstep(absTol, relTol);
    typedef boost::numeric::odeint::runge_kutta_dopri5<state_typeGPU, value_type, state_typeGPU, value_type> stepperType;

    // Define the force function
    abcSystem abc(extent);
    // lorenz_system lorenz(extent);

    // Calculate using an adaptive stepper - do smaller step chunks to allow for error printing

    double t = initTime;
    double tDelta = (finalTime - initTime)/1000.0;

    std::cout << "Transferred - computing...";
    start = std::chrono::high_resolution_clock::now();

    while(t < finalTime)
    {

        // boost::numeric::odeint::integrate_adaptive(stepstep, abc, dX, t, t + tDelta, 0.1);
        boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(1e-010, 1e-010, stepperType()), abc, std::make_pair(dX.begin(), dX.begin() + 3 * extent), t, t+tDelta, 0.0001);

        t += tDelta;

        std::cout << "Finished " << t/(finalTime) * 100. << "% of the integration." << std::endl;

    }

    end = std::chrono::high_resolution_clock::now();
    std::cout << "Complete. Time required: " <<  std::chrono::duration<double, std::milli>(end - start).count()/1000. << "s." << std::endl;

    // Copy from device to host
    std::cout << "Copying data from device to host...";
    start = std::chrono::high_resolution_clock::now();
    x = dX;
    end = std::chrono::high_resolution_clock::now();
    std::cout << "done. Time required: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms." << std::endl;

    // Now re-assign in reverse
    std::cout << "Transferring data into CPU array...";
    start = std::chrono::high_resolution_clock::now();
    // #pragma omp parallel for
    for (int i = 0; i < (int)xSize-1; i++)
    {  
        for (int j = 0; j < (int)ySize-1; j++)
        {
            for (int k = 0; k < (int)zSize-1; k++)
            {

                xyz tmp;

                tmp.x = x[i * xSize * ySize + j * ySize + k];
                tmp.y = x[extent + i * xSize * ySize + j * ySize + k];
                tmp.z = x[2*extent + i * xSize * ySize + j * ySize + k];

                finalPos[i][j][k] = tmp;

            } // k

        } // j
        
    } // i
    end = std::chrono::high_resolution_clock::now();
    std::cout << "done. Time required: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms." << std::endl;
    xyz temp = finalPos[3][3][3];
    std::cout << "here" << std::endl; std::cout.flush();
    std::cout << "Random position GPU: " << temp.x << " " << temp.y << " " << temp.z << std::endl; std::cout.flush();

}