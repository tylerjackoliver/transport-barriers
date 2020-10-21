#ifndef __POSITION_H__
#define __POSITION_H__

#include "Field.hpp"

    /* @brief Field of particle positions
    * @tparam T numeric type of the field
    */

namespace LCS
{

    template <class Type>
    class position : public field<Type>
    {

        private:
            typedef boost::multi_array_types::extent_range range;
            double xStep = 0.0;
            double yStep = 0.0;
            double zStep = 0.0;
            double xMin = 0.0;
            double xMax = 0.0;
            double yMin = 0.0;
            double yMax = 0.0;
            double zMin = 0.0;
            double zMax = 0.0;

            template <typename T> int sgn(T val)
            {
                return (val > 0) - (val < 0);              
            }

            void integrate(std::vector<Type>& x)
            {
                typedef std::vector<Type> state_type;
                int sign = sgn(this->getFinalTime() - this->getInitialTime());
                boost::numeric::odeint::bulirsch_stoer<state_type> bulirsch_stepper(this->getAbsTol(), this->getRelTol());
                boost::numeric::odeint::integrate_adaptive(bulirsch_stepper, dynSystem, x, this->getInitialTime(), this->getFinalTime(), sign*.01, abcFlowObserver);
            }

            /* @brief Computes the eigenvectors and eigenvalues for the flow using finite-differencing.
            */
            void _computeFlowFiniteDifferencing()
            {

                unsigned prog = 0;

                #pragma omp parallel for shared(prog)
                for (unsigned i = 0; i < this->nx_; ++i)
                {
                    Point<Type> xNext, xPrev, yNext, yPrev, zNext, zPrev;
                    Point<Type> x0Next, x0Prev, y0Next, y0Prev, z0Next, z0Prev;
                    Eigen::Matrix<Type, 3, 3> deformation, cauchy_green;
                    for (unsigned j = 0; j < this->ny_; ++j)
                    {
                        for (unsigned k = 0; k < this->nz_; ++k){
                            
                            xPrev = this->getUpdatedValue(i-1, j, k);
                            xNext = this->getUpdatedValue(i+1, j, k);
                            
                            yPrev = this->getUpdatedValue(i,j-1,k);
                            yNext = this->getUpdatedValue(i, j+1, k);
                            
                            zPrev = this->getUpdatedValue(i, j, k-1);
                            zNext = this->getUpdatedValue(i, j, k+1);
                            
                            x0Prev = this->getValue(i-1, j, k);
                            x0Next = this->getValue(i+1, j, k);
                            
                            y0Prev = this->getValue(i, j-1, k);
                            y0Next = this->getValue(i, j+1, k);
                            
                            z0Prev = this->getValue(i, j, k-1);
                            z0Next = this->getValue(i, j, k+1);

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
                            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(cauchy_green);
                            this->eigenvalueField_[i][j][k] = solver.eigenvalues().maxCoeff();
                            this->eigenvectorField_[i][j][k] = solver.eigenvectors().col(2).real();

                        }
                        
                    }
                
                    prog++;
                    std::cout << "Completed processing " << double(prog)/this->nx_ * 100. << "% of the eigen field." << std::endl;
                }
            }

            void _computeFlowAuxiliaryGrid()
            {
                unsigned prog = 0;

                // Get auxiliary grid sizing
                double auxGridSizeX = this->xStep * this->auxiliaryGridSizingFactor;
                double auxGridSizeY = this->yStep * this->auxiliaryGridSizingFactor;
                double auxGridSizeZ = this->zStep * this->auxiliaryGridSizingFactor;

                #pragma omp parallel for shared(prog)
                for (int i = 0; i < this->nx_+2; ++i) //  Rebase to zero due to OpenMP not supporting -ve sentinels
                {
                    for (int j = 0; j < this->ny_+2; ++j)
                    {
                        for (int k = 0; k < this->nz_+2; ++k){
                            Point<Type> xNext, xPrev, yNext, yPrev, zNext, zPrev;
                            Point<Type> x0Next, x0Prev, y0Next, y0Prev, z0Next, z0Prev;
                            Point<Type> referencePoint;
                            Eigen::Matrix<Type, 3, 3> deformation, cauchy_green;

                            std::vector<double> reference(3), xPlus(3), xMinus(3), yPlus(3), yMinus(3), zPlus(3), zMinus(3);
                            referencePoint = this->getValue(i-1, j-1, k-1);

                            reference[0] = referencePoint.x; reference[1] = referencePoint.y;
                            reference[2] = referencePoint.z;

                            xPlus = reference; xPlus[0] += auxGridSizeX; x0Next.vectorToPoint(xPlus);
                            xMinus = reference; xMinus[0] -= auxGridSizeX; x0Prev.vectorToPoint(xMinus);
                            yPlus = reference; yPlus[1] += auxGridSizeY; y0Next.vectorToPoint(yPlus);
                            yMinus = reference; yMinus[1] -= auxGridSizeY; y0Prev.vectorToPoint(yMinus);
                            zPlus = reference; zPlus[2] += auxGridSizeZ; z0Next.vectorToPoint(zPlus);
                            zMinus = reference; zMinus[2] -= auxGridSizeZ; z0Prev.vectorToPoint(zMinus);

                            this->integrate(xPlus); xNext.vectorToPoint(xPlus);  // Explicit casting std::vector -> struct Point
                            this->integrate(xMinus); xPrev.vectorToPoint(xMinus);
                            this->integrate(yPlus); yNext.vectorToPoint(yPlus);
                            this->integrate(yMinus); yPrev.vectorToPoint(yMinus);
                            this->integrate(zPlus); zNext.vectorToPoint(zPlus);
                            this->integrate(zMinus); zPrev.vectorToPoint(zMinus);

                            // deformation tensor 
                            deformation(0,0) = (xNext.x-xPrev.x) / (2. * auxGridSizeX);
                            deformation(0,1) = (yNext.x-yPrev.x) / (2. * auxGridSizeY);
                            deformation(0,2) = (zNext.x-zPrev.x) / (2. * auxGridSizeZ);
                            
                            deformation(1,0) = (xNext.y-xPrev.y) / (2. * auxGridSizeX);
                            deformation(1,1) = (yNext.y-yPrev.y) / (2. * auxGridSizeY);
                            deformation(1,2) = (zNext.y-zPrev.y) / (2. * auxGridSizeZ);

                            deformation(2,0) = (xNext.z-xPrev.z) / (2. * auxGridSizeX);
                            deformation(2,1) = (yNext.z-yPrev.z) / (2. * auxGridSizeY);
                            deformation(2,2) = (zNext.z-zPrev.z) / (2. * auxGridSizeZ);
                            
                            cauchy_green = deformation.transpose() * deformation;
                            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(cauchy_green);
                            this->eigenvalueField_[i-1][j-1][k-1] = solver.eigenvalues().maxCoeff();
                            this->eigenvectorField_[i-1][j-1][k-1] = solver.eigenvectors().col(2).real();
                        }
                        
                    }
                
                    prog++;
                    std::cout << "Completed processing " << double(prog)/this->nx_ * 100. << "% of the eigen field." << std::endl;
                }
            }

        public:

            boost::multi_array<Point<Type>, 3> updatedPosition_;
            boost::multi_array<double, 3> eigenvalueField_;
            boost::multi_array<Eigen::Vector3d, 3> eigenvectorField_;
            double initTime_;
            double finalTime_;
            double absTol_ = 1e-08; // Integrator settings
            double relTol_ = 1e-08; // Integrator settings
            double auxiliaryGridSizingFactor = 0.25;
            bool useAuxiliaryGrid = true;

            /* Constructor for initialising the position field
            * @param[in] nx The number of points in the \f$x\f$-direction
            * @param[in] ny The number of points in the \f$y\f$-direction
            * @param[in] nz The number of points in the \f$z\f$-direction
            */

            position(unsigned nx, unsigned ny, unsigned nz): field<Type>(nx, ny, nz)
            {
                typedef boost::multi_array<Point<Type>, 3> tempPoint;
                typedef boost::multi_array<double, 3> tempDouble;
                typedef boost::multi_array<Eigen::Vector3d, 3> tempVector;
                typename tempPoint::extent_gen extentsPoint;
                typename tempDouble::extent_gen extentsDouble;
                typename tempVector::extent_gen extentsVector;
                this->updatedPosition_.resize(extentsPoint[nx+2][ny+2][nz+2]);
                this->eigenvalueField_.resize(extentsDouble[nx+2][ny+2][nz+2]);
                this->eigenvectorField_.resize(extentsVector[nx+2][ny+2][nz+2]);
                this->updatedPosition_.reindex(-1);
                this->eigenvalueField_.reindex(-1);
                this->eigenvectorField_.reindex(-1);
            }

            /* Get the position of an entry in the field.

            @param[in] i Index of the grid in the x-direction
            @param[in] j Index of the grid in the y-direction
            @param[in] k Index of the grid in the z-direction

            @returns Point<Type> Value of the grid at that point
            */
            Point<Type> getValue(int i, int j, int k) // int as may want to access borders
            {
                return this->data_[i][j][k];
            }

            Type getXMax()
            {
                return this->xMax;
            }

            Type getXMin()
            {
                return this->xMin;
            }

            Type getYMin()
            {
                return this->yMin;
            }

            Type getYMax()
            {
                return this->yMax;
            }

            Type getZMin()
            {
                return this->zMin;
            }

            Type getZMax()
            {
                return this->zMax;
            }

            /* Get updated position of an entry in the field.

            @param[in] i Index of the grid in the x-direction
            @param[in] j Index of the grid in the y-direction
            @param[in] k Index of the grid in the z-direction

            @returns Point<Type> Value of the advected grid at that point
            */
            Point<Type> getUpdatedValue(int i, int j, int k) // int as may want to access borders
            {
                return this->updatedPosition_[i][j][k];
            }

            /* @brief Set data values of the field using \f$x\f$, \f$y\f$, \f$z\f$ ranges
            *
            * @param[in] xrange Vector that contains \f$x\f$ coordinates of the grid
            * @param[in] yrange Vector that contains \f$y\f$ coordinates of the grid
            * @param[in] zrange Vector that contains \f$z\f$ coordinates of the grid
            */

            void setAll(const std::vector<Type>& xrange, const std::vector<Type>& yrange, const std::vector<Type>& zrange)
            {
                // Make sure the sizes match - note the indexing is -1
                if (xrange.size()!=this->nx_+2 || yrange.size()!=this->ny_+2 \
                        || zrange.size()!=this->nz_+2)
                {
                    throw std::domain_error("Input ranges in setAll do not match.");
                }

                this->xMin = xrange[0];
                this->xMax = xrange[xrange.size()-1];

                this->yMin = yrange[0];
                this->yMax = yrange[yrange.size()-1];

                this->zMin = zrange[0];
                this->zMax = zrange[zrange.size()-1];

                // Add data to array
                #pragma omp parallel for
                for (int i = -1; i < (int)this->nx_+1; ++i)
                {   for(int j = -1; j < (int)this->ny_+1; ++j)
                    {   for (int k = -1; k < (int)this->nz_+1; ++k)
                        {
                            this->data_[i][j][k].x = xrange[i+1];
                            this->data_[i][j][k].y = yrange[j+1];
                            this->data_[i][j][k].z = zrange[k+1];
                        }
                    }
                }
            }

            /* @brief Set data values of the field using \f$x\f$, \f$y\f$, \f$z\f$ coordinates
            *
            * @param[in] xmin Minimum \f$x\f$ coordinate
            * @param[in] xmax Maximum \f$x\f$ coordinate
            * @param[in] ymin Minimum \f$y\f$ coordinate
            * @param[in] ymax Maximum \f$y\f$ coordinate
            * @param[in] zmin Minimum \f$z\f$ coordinate
            * @param[in] zmax Maximum \f$z\f$ coordinate
            */
            void setAll(Type xmin, Type xmax, Type ymin, Type ymax, Type zmin, Type zmax)
            {
                std::vector<Type> xrange(this->nx_+2, 0.), yrange(this->ny_+2, 0.), zrange(this->nz_+2, 0.);

                int i = -1;
                int j = -1;
                int k = -1; // due to -1 extent

                // Fill in values with uniform step
                if (this->nx_ < 2)
                {
                    Type gap = (xmax - xmin) / (this->nx_+1);
                    this->xStep = gap;
                    i = 0;
                    std::generate(xrange.begin(), xrange.end(), \
                        [&]{ return xmin + (i++) * gap;});
                } else
                {
                    this->xStep = (xmax-xmin) / (this->nx_-1);
                    std::generate(xrange.begin(), xrange.end(), \
                        [&]{ return xmin + (i++) * (xmax-xmin) / (this->nx_-1);});
                }

                if (this->ny_ < 2)
                {
                    Type gap = (ymax - ymin) / (this->ny_+1);
                    this->yStep = gap;
                    j = 0;
                    std::generate(yrange.begin(), yrange.end(), 
                        [&]{ return ymin + (j++) * gap;});
                } else
                {
                    this->yStep = (ymax-ymin) / (this->ny_-1);
                    std::generate(yrange.begin(), yrange.end(), \
                        [&]{ return ymin + (j++) * (ymax-ymin) / (this->ny_-1);});
                }

                // Correction for if we're on a single z-plane (planar flow)
                if (this->nz_ < 2)
                {
                    Type gap = (zmax - zmin) / (this->nz_+1);
                    this->zStep = gap;
                    k = 0;
                    std::generate(zrange.begin(), zrange.end(), 
                        [&]{ return zmin + (k++) * gap;});
                } else
                {
                    this->zStep = (zmax-zmin) / (this->nz_-1);
                    std::generate(zrange.begin(), zrange.end(), \
                        [&]{ return zmin + (k++) * (zmax-zmin) / (this->nz_-1);});
                }

                setAll(xrange, yrange, zrange);
            }

            /* @brief Set the initial time associated with the particles of the field
            *
            * @param[in] initTime Initial time of the field
            */
            void setInitialTime(double initTime)
            {
                this->initTime_ = initTime;
            }

            /* @brief Set the final time associated with the advection of the particles of the field.
            *
            * @param[in] finalTime Final time associated with the advection of the particles of the field.
            */
            void setFinalTime(double finalTime)
            {
                this->finalTime_ = finalTime;
            }

            /* @brief Get the initial time associated with the particles of the field.
            *
            * @returns The initial time associated with the particles of the field.
            */
            double getInitialTime()
            {
                return this->initTime_;
            }

            /* @brief Get the final time associated with the particles of the field.
            *
            * @returns The final time associated with the particles of the field.
            */
            double getFinalTime()
            {
                return this->finalTime_;
            }

            /* @brief Get the absolute integration tolerance for advecting the field.
            *
            * @returns The absolute integration tolerance for advecting the field.
            */
            double getAbsTol()
            {
                return this->absTol_;
            }

            /* @brief Get the relative integration tolerance for advecting the field.
            *
            * @returns The relative integration tolerance for advecting the field.
            */
            double getRelTol()
            {
                return this->relTol_;
            }

            /* @brief Set the absolute integration tolerance for advecting the field.
            *
            * @param[in] absTol The absolute integration tolerance for advecting the field.
            */
            void setAbsTol(double absTol)
            {
                this->absTol_ = absTol;
            }

            /* @brief Set the relative integration tolerance for advecting the field.
            *
            * @param[in] relTol The relative integration tolerance for advecting the field.
            */
            void setRelTol(double relTol)
            {
                this->relTol_ = relTol;
            }

            double getXStep() const
            {
                return this->xStep;
            }

            double getYStep() const
            {
                return this->yStep;
            }

            double getZStep() const
            {
                return this->zStep;
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
                            typedef std::vector<Type> state_type;
                            state_type x(3);
                            x[0] = this->data_[i][j][k].x;
                            x[1] = this->data_[i][j][k].y;
                            x[2] = this->data_[i][j][k].z;

                            // Get the sign of the time period for the numerical integration
                            int sign = this->sgn(this->finalTime_ - this->initTime_);

                            // Advect forward

                            boost::numeric::odeint::bulirsch_stoer<state_type> bulirsch_stepper(absTol, relTol);
                            boost::numeric::odeint::integrate_adaptive(bulirsch_stepper, dynSystem, x, this->initTime_, this->finalTime_, sign*.01, abcFlowObserver);
                            
                            // x now holds the final state; use this to update
                            
                            this->updatedPosition_[i][j][k].x = x[0];
                            this->updatedPosition_[i][j][k].y = x[1];
                            this->updatedPosition_[i][j][k].z = x[2];
                        } // k
                    } // j
                    #pragma omp atomic
                    prog++;
                    std::cout << "Completed advecting " << double(prog)/(this->nx_+1) * 100. << "% of the flow." << std::endl;
                } // i
            } // void advectPosition()

            void advectPositionGPU(double absTol, double relTol)
            {
                advectPositionGPUDriver(this->data_, updatedPosition_, initTime_, finalTime_, absTol, relTol);
                Type temp = updatedPosition_[3][3][3];
            }

            /* @brief Get the number of points in the 'x'-direction
            *
            * @returns The number of points in the 'x'-direction
            */
            auto getXExtent() const
            {
                return this->nx_;
            }

            /* @brief Get the number of points in the 'y'-direction
            *
            * @returns The number of points in the 'y'-direction
            */
            auto getYExtent() const
            {
                return this->ny_;
            }

            /* @brief Get the number of points in the 'z'-direction
            *
            * @returns The number of points in the 'z'-direction
            */
            auto getZExtent() const
            {
                return this->nz_;
            }

            /* @brief Write the updated position to a file "field.data"
            */
            void writeToFile() const
            {
                std::ofstream output;
                output.open("field.data", std::ofstream::out);
                
                for(unsigned i = 0; i < this->nx_; ++i)
                {
                    for (unsigned j=0; j<this->ny_; ++j)
                    {
                        for (unsigned k=0; k<this->nz_; ++k)
                        {
                            output << this->updatedPosition_[i][j][k].x << std::endl;
                            output << this->updatedPosition_[i][j][k].y << std::endl;
                            output << this->updatedPosition_[i][j][k].z << std::endl;
                        }
                    }
                }
                output.close();
            }

            /* @brief Write the updated position to a given file fname.
            @param[in] fname std::string of filename to write to.
            */
            void writeToFile(std::string fName) const
            {
                std::ofstream output;
                output.open(fName);
                
                for(unsigned i = 0; i < this->nx_; ++i)
                {
                    for (unsigned j=0; j<this->ny_; ++j)
                    {
                        for (unsigned k=0; k<this->nz_; ++k)
                        {
                            output << this->updatedPosition_[i][j][k].x << std::endl;
                            output << this->updatedPosition_[i][j][k].y << std::endl;
                            output << this->updatedPosition_[i][j][k].z << std::endl;
                        }
                    }
                }
                output.close();
            }

            /* @brief Read the upadted positions from a file field.data
            */
            void readFromFile()
            {
                std::ifstream input;
                input.open("field.data");
                
                for(unsigned i = 0; i < this->nx_; ++i)
                {
                    for (unsigned j=0; j<this->ny_; ++j)
                    {
                        for (unsigned k=0; k<this->nz_; ++k)
                        {
                            input >> this->updatedPosition_[i][j][k].x;
                            input >> this->updatedPosition_[i][j][k].y;
                            input >> this->updatedPosition_[i][j][k].z;
                        }
                    }
                }
                input.close();
            }

            /* @brief Read the updated positions from a file fname
            * @param[in] fname std::string giving the name of the file to load positions from.
            */
            void readFromFile(std::string fname)
            {
                std::ifstream input;
                input.open(fname);
                
                for(unsigned i = 0; i < this->nx_; ++i)
                {
                    for (unsigned j=0; j<this->ny_; ++j)
                    {
                        for (unsigned k=0; k<this->nz_; ++k)
                        {
                            input >> this->updatedPosition_[i][j][k].x;
                            input >> this->updatedPosition_[i][j][k].y;
                            input >> this->updatedPosition_[i][j][k].z;
                        }
                    }
                }
                input.close();
            }

            /* @brief Computes the eigenvectors and eigenvalues for the flow.
            */
            void computeFlowProperties()
            {
                if (useAuxiliaryGrid) _computeFlowAuxiliaryGrid();
                else _computeFlowFiniteDifferencing();
            }

            /* @brief Getter function for the eigenvalue field
               @param[in] i x-coordinate of the desired point
               @param[in] j y-coordinate of the desired point
               @param[in] k z-coordinate of the desired point
               @returns The eigenvalue at that point
            */
           double getEigenvalue(int i, int j, int k)
           {
               return this->eigenvalueField_[i][j][k];
           }

            /* @brief Getter function for the eigenvector field
               @param[in] i x-coordinate of the desired point
               @param[in] j y-coordinate of the desired point
               @param[in] k z-coordinate of the desired point
               @returns The eigenvector at that point
            */
           Eigen::Vector3d getEigenvector(int i, int j, int k)
           {
               return this->eigenvectorField_[i][j][k];
           }
    }; // class
}; // namespace
#endif