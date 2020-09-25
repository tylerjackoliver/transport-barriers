#pragma once

#include <omp.h>
#include <boost/multi_array.hpp>
#include <eigen3/Eigen/Eigenvalues>
#include "commonGPUCPU.hpp"
#include "xyz.hpp" // struct
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>
#include <fstream>
#include "forceFunction.hpp"
#include <iostream>
#include <cmath>
#include <utility>
#include <chrono>
#include <cuda.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

// CUDA definitions
typedef double value_type;
typedef thrust::host_vector<value_type> state_typeCPU;
typedef thrust::device_vector<value_type> state_typeGPU;

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


struct planarDoubleGyre{ //  Double gyre with zero z-accelerations (so permanently confined to plane)

    struct doubleGyreFunctor
    {

        template<class T>
        __host__ __device__ 
        void operator()(T t) const{

            // Unpack the parameters
            value_type x = thrust::get<0>(t);
            value_type y = thrust::get<1>(t);
            value_type z = thrust::get<2>(t);

            value_type pi_ = std::atan(1.0) * 4.0; // 2pi

            // thrust::get<3>(t) = sqrt(3.0) * sin(z) + 1.0 * cos(y);
            // thrust::get<4>(t) = sqrt(2.0) * sin(x) + sqrt(3.0) * cos(z);
            // thrust::get<5>(t) = 1.0 * sin(y) + sqrt(2) * cos(x);

            thrust::get<3>(t) = -pi_ * std::sin(pi_*x) * std::cos(pi_*y);
            thrust::get<4>(t) = pi_ * std::cos(pi_ * x) * std::sin(pi_ * y);
            thrust::get<5>(t) = 0.0;

        }

    };

    planarDoubleGyre(size_t N) : dim_N(N){}

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

            ), doubleGyreFunctor() );

    };

    size_t dim_N;

};

void advectPositionGPUDriver(field_& initPos, field_& updatedPos, double initTime, double finalTime, \
double absTol, double relTol);

namespace LCS
{
    class field
    {

        public:

            /* @brief Constructor for the field class.
             *
             * @param nx The number of points in x
             * @param ny The number of points in y
             * @param nz The number of points in z
             */
            field(unsigned nx, unsigned ny, unsigned nz) : nx_(nx), ny_(ny), nz_(nz), data_(boost::extents[range(-1,(int)(nx_+1))][range(-1, (int)(ny_+1))][range(-1, (int)(nz_+1))])
            {

            };

            /* @brief Retrieve elements of the domain
             *
             * @param i The index of the point in x
             * @param j The index of the point in y
             * @param k The index of the point in z
             */

            xyz getValue(unsigned i, unsigned j, unsigned k)
            {
    
                return data_[i][j][k];

            }

            /* @brief Set a value of the underlying array
             *
             * @param value The value to set
             * @param i The index of the point in x
             * @param j The index of the point in y
             * @param k The index of the point in z
             */

            void setValue(xyz value, unsigned i, unsigned j, unsigned k)
            {

                data_[i][j][k] = value;

            }
            
        protected:

            typedef std::vector<double> state_type;
            unsigned nx_;
            unsigned ny_;
            unsigned nz_;
            field_ data_;

    };

    /* @brief Field of particle positions
    * @tparam T numeric type of the field
    */

    class position : public field
    {

        public:

            /* Constructor for initialising the position field
            * @param nx The number of points in the \f$x\f$-direction
            * @param ny The number of points in the \f$y\f$-direction
            * @param nz The number of points in the \f$z\f$-direction
            */

            position(unsigned nx, unsigned ny, unsigned nz): field(nx, ny, nz), updatedPosition_(extents[range(-1, (int)(this->nx_+1))][range(-1, (int)(this->ny_+1))][range(-1, (int)(this->nz_+1))]){
            
                // field_ updatedPosition_(extents[range(-1, (int)(this->nx_+1))][range(-1, (int)(this->ny_+1))][range(-1, (int)(this->nz_+1))]);

            }

            /* Get the position of an entry in the field.

            @param i Index of the grid in the x-direction
            @param j Index of the grid in the y-direction
            @param k Index of the grid in the z-direction
            */
            xyz getValue(int i, int j, int k) // int as may want to access borders
            {

                return this->data_[i][j][k]; // Due to -1 indexing

            }

            /* Get updated position of an entry in the field.

            @param i Index of the grid in the x-direction
            @param j Index of the grid in the y-direction
            @param k Index of the grid in the z-direction
            */
            xyz getUpdatedValue(int i, int j, int k) // int as may want to access borders
            {

                return this->updatedPosition_[i][j][k]; // Due to -1 indexing

            }

            /* @brief Set data values of the field using \f$x\f$, \f$y\f$, \f$z\f$ ranges
            *
            * @param xrange Vector that contains \f$x\f$ coordinates of the grid
            * @param yrange Vector that contains \f$y\f$ coordinates of the grid
            * @param zrange Vector that contains \f$z\f$ coordinates of the grid
            */

            void setAll(const std::vector<double>& xrange, const std::vector<double>& yrange, const std::vector<double>& zrange)
            {

                // Make sure the sizes match - note the indexing is -1
                if (xrange.size()!=this->nx_+2 || yrange.size()!=this->ny_+2 \
                        || zrange.size()!=this->nz_+2)
                {

                    throw std::domain_error("Input ranges in setAll do not match.");

                }

                // Add data to array
                #pragma omp parallel for
                for (int i = -1; i < (int)this->nx_+1; ++i)
                {   for(int j = -1; j < (int)this->ny_+1; ++j)
                    {   for (int k = -1; k < (int)this->nz_+1; ++k)
                        {
                            data_[i][j][k].x = xrange[i+1];
                            data_[i][j][k].y = yrange[j+1];
                            data_[i][j][k].z = zrange[k+1];
                        }
                    }
                }

            }

            /* @brief Set data values of the field using \f$x\f$, \f$y\f$, \f$z\f$ coordinates
            *
            * @param xmin Minimum \f$x\f$ coordinate
            * @param xmax Maximum \f$x\f$ coordinate
            * @param ymin Minimum \f$y\f$ coordinate
            * @param ymax Maximum \f$y\f$ coordinate
            * @param zmin Minimum \f$z\f$ coordinate
            * @param zmax Maximum \f$z\f$ coordinate
            */

            void setAll(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
            {

                std::vector<double> xrange(this->nx_+2, 0.), yrange(this->ny_+2, 0.), zrange(this->nz_+2, 0.);

                int i = -1;
                int j = -1;
                int k = -1; // due to -1 extent

                // Fill in values with uniform step
                
                std::generate(xrange.begin(), xrange.end(), \
                        [&]{ return xmin + (i++) * (xmax-xmin) / (this->nx_-1);});
                std::generate(yrange.begin(), yrange.end(), \
                        [&]{ return ymin + (j++) * (ymax-ymin) / (this->ny_-1);});
                std::generate(zrange.begin(), zrange.end(), \
                        [&]{ return zmin + (k++) * (zmax-zmin) / (this->nz_-1);});

                setAll(xrange, yrange, zrange);
                    
            }

            void setInitialTime(double initTime){

                this->initTime_ = initTime;

            }

            void setFinalTime(double finalTime){

                this->finalTime_ = finalTime;

            }

            double getInitialTime()
            {

                return this->initTime_;

            }

            double getFinalTime()
            {

                return this->finalTime_;

            }

            /* @brief Advect the initial flow coordinates forward using a Burlisch-Stoer numerical integration scheme
            *
            * @param absTol Absolute tolerance to use in the integration
            * @param relTol Relative tolerance to use in the integration
            */

            void advectPosition(double absTol, double relTol)
            {

                unsigned prog = 0;

                #pragma omp parallel for shared(prog)
                for (int i = -1; i < (int)this->nx_+1; i++)
                {   for (int j = -1; j < (int)this->ny_+1; j++)
                    {
                        for (int k = -1; k < (int)this->nz_+1; k++)
                        {

                            state_type x(3);
                            x[0] = data_[i][j][k].x;
                            x[1] = data_[i][j][k].y;
                            x[2] = data_[i][j][k].z;

                            // Advect forward
                            
                            boost::numeric::odeint::bulirsch_stoer<state_type> bulirsch_stepper(absTol, relTol);
                            boost::numeric::odeint::integrate_adaptive(bulirsch_stepper, abcFlow, x, this->initTime_, this->finalTime_, .01, abcFlowObserver);

                            // x now holds the final state; use this to update
                            
                            this->updatedPosition_[i][j][k].x = x[0];
                            this->updatedPosition_[i][j][k].y = x[1];
                            this->updatedPosition_[i][j][k].z = x[2];

                        } // k

                    } // j

                    prog++;
                    std::cout << "Completed " << double(prog)/(this->nx_+1) * 100. << "%." << std::endl;

                } // i

            } // void advectPosition()

            void advectPositionGPU(double absTol, double relTol){

                int driverVersion, runtimeVersion;
                cudaDriverGetVersion(&driverVersion); cudaRuntimeGetVersion(&runtimeVersion);
            
                std::cout << "Computing on CUDA version \t" << driverVersion << "\t and runtime version \t" << runtimeVersion << std::endl;
            
                // Get problem dimensionality
                unsigned int xSize = this->data_.shape()[0];
                unsigned int ySize = this->data_.shape()[1];
                unsigned int zSize = this->data_.shape()[2];
                // Get the size of the flattened array
                unsigned int extent = xSize * ySize * zSize;
                // Number of dimensions for the array
                unsigned int dim = 3;
                // Assign flattened array
                state_typeCPU x(dim*extent);
            
                /* Fill the device vector with the correct values.
                *  For speed, initialise _first_ on the host device (x), then transfer to the GPU in one chunk (dx)
                */
            
                std::cout << "Transferring data into CPU array..."; std::cout.flush();
                auto start = std::chrono::high_resolution_clock::now();
            
                for (int i = 0; i < (int)this->nx_+1; i++)
                {  
                    for (int j = 0; j < (int)this->ny_+1; j++)
                    {
                        for (int k = 0; k < (int)this->nz_+1; k++)
                        {
            
                            xyz tmp = this->data_[i][j][k];
            
                            x[i * zSize * ySize + j * zSize + k] = tmp.x;
                            x[extent + i * zSize * ySize + j * zSize + k] = tmp.y;
                            x[2*extent + i * zSize * ySize + j * zSize + k] = tmp.z; 
            
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
                planarDoubleGyre abc(extent);
                // lorenz_system lorenz(extent);
            
                // Calculate using an adaptive stepper - do smaller step chunks to allow for error printing
            
                double t = this->initTime_;
                double tDelta = (this->finalTime_ - this->initTime_)/1000.0;
            
                std::cout << "Transferred - computing...";
                start = std::chrono::high_resolution_clock::now();
            
                while(t < this->finalTime_)
                {
            
                    // boost::numeric::odeint::integrate_adaptive(stepstep, abc, dX, t, t + tDelta, 0.1);
                    boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(1e-010, 1e-010, stepperType()), abc, std::make_pair(dX.begin(), dX.begin() + 3 * extent), t, t+tDelta, 0.0001);
            
                    t += tDelta;

                    std::cout << "Finished " << t/(this->finalTime_) * 100. << "% of the integration." << std::endl;
            
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
                for (int i = 0; i < this->nx_+1; i++)
                {  
                    for (int j = 0; j < this->ny_+1; j++)
                    {
                        for (int k = 0; k < this->nz_+1; k++)
                        {
            
                            xyz tmp;
                            xyz zTemp = data_[i][j][k];
            
                            tmp.x = x[i * zSize * ySize + j * zSize + k];
                            tmp.y = x[extent + i * zSize * ySize + j * zSize + k];
                            tmp.z = zTemp.z;

                            updatedPosition_[i][j][k] = tmp;
                            std::cout << tmp.x << " " << tmp.y << " " << tmp.z << std::endl; std::cout.flush();
            
                        } // k
            
                    } // j
                    
                } // i
            
                end = std::chrono::high_resolution_clock::now();
            
            }

            unsigned getXExtent()
            {

                return nx_;

            }

            unsigned getYExtent()
            {

                return ny_;

            }

            unsigned getZExtent()
            {

                return nz_;

            }

            xyz getOriginalPosition(int i, int j, int k)
            {

                return this->data_[i][j][k];

            }

            void writeToFile(){

                std::ofstream integrationOutput;
                integrationOutput.open("DGintegration.data2", std::ofstream::out);

                std::cout << "Writing output to file...";
                for(unsigned i = 0; i < this->nx_; ++i){

                    for(unsigned j = 0; j < this->ny_; ++j)
                    {
                        for (unsigned k = 0; k < this->nz_; k++)
                        {

                            integrationOutput << this->updatedPosition_[i][j][k].x << std::endl;
                            integrationOutput << this->updatedPosition_[i][j][k].y << std::endl;
                            integrationOutput << this->updatedPosition_[i][j][k].z << std::endl;

                        }


                    }

                }

                std::cout << "done.";
                integrationOutput.flush();
                integrationOutput.close();

            }

            void loadFromFile(const char* input)
            {

                std::ifstream fileIn;
                fileIn.open(input);

                std::cout << "Reading flow-field data from text file...";

                for(unsigned i = 0; i < this->nx_; ++i){

                    for(unsigned j = 0; j < this->ny_; ++j)
                    {
                        for (unsigned k = 0; k < this->nz_; k++)
                        {

                            fileIn >> this->updatedPosition_[i][j][k].x;
                            fileIn >> this->updatedPosition_[i][j][k].y;
                            fileIn >> this->updatedPosition_[i][j][k].z;

                        }

                    }

                }
                
                std::cout << "done";

            }


            field_ updatedPosition_;
            double initTime_;
            double finalTime_;

        protected:


    }; // class

    class FTLE : public field
    {

        public:

            FTLE(position& pos): field(pos.getXExtent(), pos.getYExtent(), pos.getZExtent()),
                pos_(pos), ftleField_(extents[range(-1, (int)(pos.getXExtent()+1))][range(-1, (int)(pos.getYExtent()+1))][range(-1, (int)(pos.getZExtent()+1))]){
                };
            
            position& pos_;
            doubleField_ ftleField_;

            void computeFTLE(){ 
                
                unsigned prog = 0;
                double coeff = 1./fabs(pos_.initTime_ - pos_.finalTime_);

                #pragma omp parallel for shared(prog)
                for (unsigned i = 0; i < this->nx_; ++i)
                {
                    xyz xNext, xPrev, yNext, yPrev, zNext, zPrev;
                    xyz x0Next, x0Prev, y0Next, y0Prev, z0Next, z0Prev;
                    Eigen::Matrix<double, 3, 3> deformation, cauchy_green;
                    for (unsigned j = 0; j < this->ny_; ++j)
                    {
                        for (unsigned k = 0; k < this->nz_; ++k){
                            
                            xPrev = pos_.getUpdatedValue(i-1, j, k);
                            xNext = pos_.getUpdatedValue(i+1, j, k);
                            
                            yPrev = pos_.getUpdatedValue(i,j-1,k);
                            yNext = pos_.getUpdatedValue(i, j+1, k);
                            
                            zPrev = pos_.getUpdatedValue(i, j, k-1);
                            zNext = pos_.getUpdatedValue(i, j, k+1);
                            
                            x0Prev = pos_.getValue(i-1, j, k);
                            x0Next = pos_.getValue(i+1, j, k);
                            
                            y0Prev = pos_.getValue(i, j-1, k);
                            y0Next = pos_.getValue(i, j+1, k);
                            
                            z0Prev = pos_.getValue(i, j, k-1);
                            z0Next = pos_.getValue(i, j, k+1);

                            // deformation tensor 
                            deformation(0,0) = (xNext.x-xPrev.x) / (x0Next.x-x0Prev.x);
                            deformation(0,1) = (yNext.x-yPrev.x) / (y0Next.y-y0Prev.y);
                            deformation(0,2) = (zNext.x-zPrev.x) / (z0Next.z-z0Prev.z);
                            
                            deformation(1,0) = (xNext.y-xPrev.y) / (x0Next.x-x0Prev.x);
                            deformation(1,1) = (yNext.y-yPrev.y) / (y0Next.y-y0Prev.y);
                            deformation(1,2) = (zNext.y-zPrev.y) / (z0Next.z-z0Prev.z);

                            deformation(2,0) = (xNext.z-xPrev.z) / (x0Next.x-x0Prev.x);
                            deformation(2,1) = (yNext.z-yPrev.z) / (y0Next.y-y0Prev.y);
                            deformation(2,2) = (zNext.z-zPrev.z) / (z0Next.z-z0Prev.z);

                            cauchy_green = deformation.transpose() * deformation;
                            auto eivals = cauchy_green.template selfadjointView<Eigen::Lower>().eigenvalues();
                            ftleField_[i][j][k] = coeff * std::log(eivals.maxCoeff());

                        }
                        
                    }
                
                    prog++;
                    if(omp_get_thread_num() == 0)
                    std::cout << "Completed processing " << double(prog)/this->nx_ * 100. << "% of the FTLE." << std::endl;

                }

            }

            void writeFTLE()
            {

                std::ofstream output;
                output.open("DGftleField.data2");

                std::cout << "Writing FTLE values to file..."; std::cout.flush();

                for(unsigned i = 0; i < this->nx_; ++i)
                {

                    for(unsigned j = 0; j < this->ny_; ++j)
                    {

                        for(unsigned k = 0; k < this->nz_; k+=500)
                        {

                            output << ftleField_[i][j][k] << std::endl;

                        }

                    }

                }

                std::cout << "done. Time required: ";

                output.close();

            }

    };

    class helicity : public field
    {
        public:

            helicity(position& pos): field(pos.getXExtent(), pos.getYExtent(), pos.getZExtent()),
                pos_(pos), helicityField_(evExtents[range(0, (int)(pos.getXExtent()))][range(0, (int)(pos.getYExtent()))][range(0, (int)(pos.getZExtent()))]),
                dominantField_(evExtents[range(0, (int)(pos.getXExtent()))][range(0, (int)(pos.getYExtent()))][range(0, (int)(pos.getZExtent()))]){

                    evComputed = 0;

                };

            void getDominantVectorField(){

                unsigned prog = 0;
                std::cout << "Computing the dominant vector field...";

                #pragma omp parallel for shared(prog)
                for (unsigned i = 0; i < this->nx_; ++i)
                {
                    xyz xNext, xPrev, yNext, yPrev, zNext, zPrev;
                    xyz x0Next, x0Prev, y0Next, y0Prev, z0Next, z0Prev;
                    Eigen::Matrix<double, 3, 3> deformation, cauchy_green;

                    for (unsigned j = 0; j < this->ny_; ++j)
                    {
                        for (unsigned k = 0; k < this->nz_; ++k){
                            
                            xPrev = pos_.getUpdatedValue(i-1, j, k);
                            xNext = pos_.getUpdatedValue(i+1, j, k);
                            
                            yPrev = pos_.getUpdatedValue(i,j-1,k);
                            yNext = pos_.getUpdatedValue(i, j+1, k);
                            
                            zPrev = pos_.getUpdatedValue(i, j, k-1);
                            zNext = pos_.getUpdatedValue(i, j, k+1);
                            
                            x0Prev = pos_.getValue(i-1, j, k);
                            x0Next = pos_.getValue(i+1, j, k);
                            
                            y0Prev = pos_.getValue(i, j-1, k);
                            y0Next = pos_.getValue(i, j+1, k);
                            
                            z0Prev = pos_.getValue(i, j, k-1);
                            z0Next = pos_.getValue(i, j, k+1);

                            // deformation tensor 
                            deformation(0,0) = (xNext.x-xPrev.x) / (x0Next.x-x0Prev.x);
                            deformation(0,1) = (yNext.x-yPrev.x) / (y0Next.y-y0Prev.y);
                            deformation(0,2) = (zNext.x-zPrev.x) / (z0Next.z-z0Prev.z);
                            
                            deformation(1,0) = (xNext.y-xPrev.y) / (x0Next.x-x0Prev.x);
                            deformation(1,1) = (yNext.y-yPrev.y) / (y0Next.y-y0Prev.y);
                            deformation(1,2) = (zNext.y-zPrev.y) / (z0Next.z-z0Prev.z);

                            deformation(2,0) = (xNext.z-xPrev.z) / (x0Next.x-x0Prev.x);
                            deformation(2,1) = (yNext.z-yPrev.z) / (y0Next.y-y0Prev.y);
                            deformation(2,2) = (zNext.z-zPrev.z) / (z0Next.z-z0Prev.z);

                            cauchy_green = deformation.transpose() * deformation;

                            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> es; // CGST is positive definite + symmetric
                            es.compute(cauchy_green);

                            // SelfAdjoint sorts the eigenvectors in ascending order, so the dominant eigenvalue is the last one!

                            dominantField_[i][j][k] = es.eigenvectors().col(2); // Should be three eigenvectors!

                        }
                        
                    }
                
                    prog++;
                    if(omp_get_thread_num() == 0)
                    std::cout << "Completed processing " << double(prog)/this->nx_ * 100. << "% of the eigenvector field." << std::endl;

                    evComputed = 1;

                }

            }

            void computeHelicity(){

                if (!evComputed){

                    std::cout << "Dominant eigenvector field has not yet been computed. Computing now..." << std::endl;
                    getDominantVectorField();

                }

                // H(eta) = curl(eta).dot(eta)

                // Go through the field, compute the curl (essentially finite-differencing the vector field), and then 
                // dot it with the other

                std::cout << "Computing the helicity of the velocity field..." << std::endl;

                unsigned prog = 0;

                // #pragma omp parallel for shared(prog)
                for (unsigned i = 0; i < this->nx_; ++i){

                    Eigen::Vector3d eta(3), curl(3), xNext, xPrev, yNext, yPrev, zNext, zPrev;
                    xyz x0Next, x0Prev, y0Next, y0Prev, z0Next, z0Prev;

                    for (unsigned j = 0; j < this->ny_; ++j)
                    {

                        for (unsigned k = 0; k < this -> nz_; ++k)
                        {

                            // Get values

                            eta = dominantField_[i][j][k];

                            xNext = i != nx_-1 ? dominantField_[i+1][j][k] : dominantField_[i][j][k];
                            x0Next = i != nx_-1 ? pos_.getUpdatedValue(i+1, j, k) : pos_.getUpdatedValue(i, j, k);

                            xPrev = i != 0 ? dominantField_[i-1][j][k] : dominantField_[i][j][k];
                            x0Prev = i != 0 ? pos_.getUpdatedValue(i-1, j, k) : pos_.getUpdatedValue(i, j, k);

                            yNext = j != ny_-1 ? dominantField_[i][j+1][k] : dominantField_[i][j][k];
                            y0Next = j != ny_-1 ? pos_.getUpdatedValue(i, j+1, k) : pos_.getUpdatedValue(i, j, k);

                            yPrev = j != 0 ? dominantField_[i][j-1][k] : dominantField_[i][j][k];
                            y0Prev = j != 0 ? pos_.getUpdatedValue(i, j-1, k) : pos_.getUpdatedValue(i, j, k);

                            zNext = k != nz_-1 ? dominantField_[i][j][k+1] : dominantField_[i][j][k];
                            z0Next = k != nz_-1 ? pos_.getUpdatedValue(i, j, k+1) : pos_.getUpdatedValue(i, j, k);

                            zPrev = k != 0 ? dominantField_[i][j][k-1] : dominantField_[i][j][k];
                            z0Prev = k != 0 ? pos_.getUpdatedValue(i, j, k-1) : pos_.getUpdatedValue(i, j, k);

                            // Compute curl elements -- see https://en.wikipedia.org/wiki/Curl_(mathematics)#Usage

                            double dfzdy = (yNext(2) - yPrev(2)) / (y0Next.y - y0Prev.y);
                            double dfydz = (zNext(1) - zPrev(1)) / (z0Next.z - z0Prev.z);

                            double dfxdz = (zNext(0) - zPrev(0)) / (z0Next.z - z0Prev.z);
                            double dfzdx = (xNext(2) - xPrev(2)) / (x0Next.x - x0Prev.x);

                            double dfydx = (xNext(1) - xPrev(1)) / (x0Next.x - x0Prev.x);
                            std::cout << "dfydx: " << dfydx << std::endl;
                            double dfxdy = (yNext(0) - yPrev(0)) / (y0Next.y - y0Prev.y);
                            std::cout << "dfxdy: " << dfzdx << std::endl;

                            curl(0) = (dfzdy - dfydz);
                            curl(1) = (dfxdz - dfzdx);
                            curl(2) = (dfydx - dfxdy);

                            helicityField_[i][j][k] = fabs(curl.dot(eta));

                        }

                    }

                    prog++;
                    if(omp_get_thread_num() == 0)
                        std::cout << "Completed processing " << (double)prog / (nx_-1) * 100.0 << "% of the helicity field." << std::endl;

                }

            }

            void writeHelicityBelowThreshold(double eps)
            {

                std::ofstream output;
                output.open("isStrainline");

                for (unsigned i = 0; i < this->nx_; ++i)
                {

                    for (unsigned j = 0; j < this->ny_; ++j)
                    {

                        for (unsigned k = 0; k < this->nz_; ++k)
                        {

                            if (helicityField_[i][j][k] < eps)
                            {

                                output << 1 << std::endl;

                            } else 
                            {

                                output << 0 << std::endl;

                            }

                        }

                    }

                }

                output.close();

            }

            void writeHelicityField()
            {

                std::ofstream output;
                output.open("helicity");

                for (unsigned i = 0; i < this->nx_; ++i)
                {

                    for (unsigned j = 0; j < this->ny_; ++j)
                    {

                        for (unsigned k = 0; k < this->nz_; ++k)
                        {

                            output << helicityField_[i][j][k] << std::endl;

                        }

                    }

                }

                output.close();


            }

            std::vector<xyz> getHelicityBelowThreshold(double eps)
            {

                std::vector<xyz> output; // Actually stores i, j, k

                for (unsigned i = 0; i < this->nx_; ++i)
                {

                    for (unsigned j = 0; j < this->ny_; ++j)
                    {

                        for (unsigned k = 0; k < this->nz_; ++k)
                        {

                            if (helicityField_[i][j][k] <= eps)
                            {
                            
                                xyz tmp;
                                tmp.x = i; tmp.y = j; tmp.z = k;

                                output.push_back(tmp);
                            
                            }

                        }

                    }

                }

            }

            std::vector<std::vector<xyz>> propagateIndices(std::vector<xyz> indices, double eps)
            {

                // Get indices of points and their neighbours
                // Get the initial conditions in x, y, z (updated grid)
                // Step-through one-by-one for every trajectory to get eta3, eta3_dot
                // Stop when H3 > eps

                std::cout << "~~~~~~~~~~~~~ Strainline propagation ~~~~~~~~~~~~~~" << std::endl;
                std::cout << std::endl;
                std::cout << "There are " << indices.size() << "candidate points." << std::endl;

                std::vector<std::vector<xyz>> output;
                
                for (unsigned seedNum = 0; seedNum < indices.size(); ++seedNum)
                {

                    unsigned i = indices.x;
                    unsigned j = indices.y;
                    unsigned k = indices.z;

                    // Get the points that corresponds to

                    xyz initPosXyz = pos_.getUpdatedValue(i, j, k);

                    // Now get its neighbours

                    xyz x0NextXyz = i != nx_-1 ? pos_.getUpdatedValue(i+1, j, k) : pos_.getUpdatedValue(i, j, k);
                    xyz x0PrevXyz = i != 0 ? pos_.getUpdatedValue(i-1, j, k) : pos_.getUpdatedValue(i, j, k);

                    xyz y0NextXyz = j != ny_-1 ? pos_.getUpdatedValue(i, j+1, k) : pos_.getUpdatedValue(i, j, k);
                    xyz y0PrevXyz = j != 0 ? pos_.getUpdatedValue(i, j-1, k) : pos_.getUpdatedValue(i, j, k);

                    xyz z0NextXyz = k != nz_-1 ? pos_.getUpdatedValue(i, j, k+1) : pos_.getUpdatedValue(i, j, k);
                    xyz z0PrevXyz = k != 0 ? pos_.getUpdatedValue(i, j, k-1) : pos_.getUpdatedValue(i, j, k);

                    // Get its initial helicity

                    double initHelicity = helicityField_[i][j][k];

                    // Print status message

                    std::cout << "Computing candidate point number " << seedNum << std::endl;

                    // Initialise integrator

                    typedef boost::numeric::odeint::runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double> stepperType;
                
                    // Initialise ICs

                    std::vector<double> initPos(12), xNext(12), xPrev(12), yNext(12), yPrev(12), zNext(12), zPrev(12);

                    initPos = {initPosXyz.x, initPosXyz.y, initPosXyz.z, 1, 0, 0, 0, 1, 0, 0, 0, 1};

                    xNext = {x0NextXyz.x, x0NextXyz.y, x0NextXyz.z, 1, 0, 0, 0, 1, 0, 0, 0, 1}; // + identity matrix
                    xPrev = {x0PrevXyz.x, x0PrevXyz.y, x0PrevXyz.z, 1, 0, 0, 0, 1, 0, 0, 0, 1};

                    yNext = {y0NextXyz.x, y0NextXyz.y, y0NextXyz.z, 1, 0, 0, 0, 1, 0, 0, 0, 1};
                    yPrev = {y0PrevXyz.x, y0PrevXyz.y, y0PrevXyz.z, 1, 0, 0, 0, 1, 0, 0, 0, 1};

                    zNext = {z0NextXyz.x, z0NextXyz.y, z0NextXyz.z, 1, 0, 0, 0, 1, 0, 0, 0, 1};
                    zPrev = {z0PrevXyz.x, z0PrevXyz.y, z0PrevXyz.z, 1, 0, 0, 0, 1, 0, 0, 0, 1};

                    // Calculate using an adaptive stepper - do smaller step chunks to allow for error printing
                
                    double t = 0.0;
                    double tDelta = 0.1;
                    unsigned numSteps = 0;

                    double helicitySum = 0.0;
                    double helicityAvg = 0.0;

                    std::vector<xyz> trajectory;
                    trajectory.push_back(initPosXyz);
                
                    while(t < 100.)
                    {
                
                        // Integrate for a small period of time

                        boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(1e-010, 1e-010, stepperType()), abcFlowVariational, initPos, t, t+tDelta, 0.0001);

                        // Neighbouring points

                        boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(1e-010, 1e-010, stepperType()), abcFlowVariational, xNext, t, t+tDelta, 0.0001);
                        boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(1e-010, 1e-010, stepperType()), abcFlowVariational, xPrev, t, t+tDelta, 0.0001);                
                        boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(1e-010, 1e-010, stepperType()), abcFlowVariational, yNext, t, t+tDelta, 0.0001);
                        boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(1e-010, 1e-010, stepperType()), abcFlowVariational, yPrev, t, t+tDelta, 0.0001);                
                        boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(1e-010, 1e-010, stepperType()), abcFlowVariational, zNext, t, t+tDelta, 0.0001);
                        boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(1e-010, 1e-010, stepperType()), abcFlowVariational, zPrev, t, t+tDelta, 0.0001);                

                        t += tDelta;
                        numSteps++;
                        xyz vectorToXyz;
                        vectorToXyz.x = initPos[0]; vectorToXyz.y = initPos[1]; vectorToXyz = initPos[2]; 
                        trajectory.push_back(vectorToXyz);

                        // Get the dominant eigenvector for every integration (from variational equations)

                        Eigen::Vector3d eig = dominantEigenvector(initPos); 

                        Eigen::Vector3d xNextEig = dominantEigenvector(xNext); 
                        Eigen::Vector3d xPrevEig = dominantEigenvector(xPrev); 
                        Eigen::Vector3d yNextEig = dominantEigenvector(yNext); 
                        Eigen::Vector3d yPrevEig = dominantEigenvector(yPrev); 
                        Eigen::Vector3d zNextEig = dominantEigenvector(zNext); 
                        Eigen::Vector3d zPrevEig = dominantEigenvector(zPrev); 

                        // Get the finite differencing in the dominant eigenvector for the curl

                        double dfzdy = (yNextEig(2) - yPrevEig(2)) / (y0NextXyz.y - y0PrevXyz.y);
                        double dfydz = (zNextEig(1) - zPrevEig(1)) / (z0NextXyz.z - z0PrevXyz.z);

                        double dfxdz = (zNextEig(0) - zPrevEig(0)) / (z0NextXyz.z - z0PrevXyz.z);
                        double dfzdx = (xNextEig(2) - xPrevEig(2)) / (x0NextXyz.x - x0PrevXyz.x);

                        double dfydx = (xNextEig(1) - xPrevEig(1)) / (x0NextXyz.x - x0PrevXyz.x);
                        double dfxdy = (yNextEig(0) - yPrevEig(0)) / (y0NextXyz.y - y0PrevXyz.y);

                        Eigen::Vector3d curl;

                        // Compute curl

                        curl(0) = (dfzdy - dfydz);
                        curl(1) = (dfxdz - dfzdx);
                        curl(2) = (dfydx - dfxdy);

                        // Compute helicity

                        double helicity = fabs(curl.dot(eig));
                        helicitySum += helicity;
                        helicityAvg = helicitySum / numSteps;

                        if (helicityAvg > eps) break;

                    }

                    std::cout << "Stopped integration at time " << t << std::endl;
                    output.push_back(trajectory);

                }

            }

            Eigen::Vector3d dominantEigenvector(std::vector<double>& input)
            {

                // Form of input is (state, state, state, derivatives->end)

                Eigen::Matrix<double, 3, 3> deformation, cgst;

                deformation(0,0) = input[3];
                deformation(0,1) = input[4];
                deformation(0,2) = input[5];

                deformation(1,0) = input[6];
                deformation(1,1) = input[7];
                deformation(1,2) = input[8];

                deformation(2,0) = input[9];
                deformation(2,1) = input[10];
                deformation(2,2) = input[11];

                cgst = deformation.transpose() * deformation;

                Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> es;
                es.compute(cgst);

                // Return the dominant eigenvector

                return es.eigenvectors.col(2);

            }

            position& pos_;
            eigenvectorField_ dominantField_;
            doubleField_ helicityField_;
            bool evComputed = 0; // Tracks whether we've computed the EV field before the helicity

    };

}; // namespace

