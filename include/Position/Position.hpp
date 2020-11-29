#ifndef __POSITION_H__
#define __POSITION_H__

#include "Field.hpp"
#include "Helpers.hpp"

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

                            double initTime = this->initTime_;
                            double finalTime = this->finalTime_;

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

                            Helpers::integrateForceFunction(xPlus, this->initTime_, this->finalTime_, this->absTol_, this->relTol_);
                            Helpers::integrateForceFunction(xMinus, this->initTime_, this->finalTime_, this->absTol_, this->relTol_);
                            Helpers::integrateForceFunction(yPlus, this->initTime_, this->finalTime_, this->absTol_, this->relTol_);
                            Helpers::integrateForceFunction(yMinus, this->initTime_, this->finalTime_, this->absTol_, this->relTol_);
                            Helpers::integrateForceFunction(zPlus, this->initTime_, this->finalTime_, this->absTol_, this->relTol_);
                            Helpers::integrateForceFunction(zMinus, this->initTime_, this->finalTime_, this->absTol_, this->relTol_);

                            xNext.vectorToPoint(xPlus);  // Explicit casting std::vector -> struct Point
                            xPrev.vectorToPoint(xMinus);
                            yNext.vectorToPoint(yPlus);
                            yPrev.vectorToPoint(yMinus);
                            zNext.vectorToPoint(zPlus);
                            zPrev.vectorToPoint(zMinus);

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
                /* Construct temporary types for the later boost::multi_arrays */
                typedef boost::multi_array<Point<Type>, 3> tempPoint;
                typedef boost::multi_array<double, 3> tempDouble;
                typedef boost::multi_array<Eigen::Vector3d, 3> tempVector;
                /* Construct temporary types for the boost extents */
                typename tempPoint::extent_gen extentsPoint;
                typename tempDouble::extent_gen extentsDouble;
                typename tempVector::extent_gen extentsVector;
                /* Resize the relevant arrays to have the same size, and same initial indices */
                this->updatedPosition_.resize(extentsPoint[nx+2][ny+2][nz+2]);
                this->eigenvalueField_.resize(extentsDouble[nx+2][ny+2][nz+2]);
                this->eigenvectorField_.resize(extentsVector[nx+2][ny+2][nz+2]);
                this->updatedPosition_.reindex(-1);
                this->eigenvalueField_.reindex(-1);
                this->eigenvectorField_.reindex(-1);
            }

            /* Function signatures */
            /* Getters/setters */
            auto getEigenvalue(int, int, int);
            Eigen::Vector3d getEigenvector(int, int, int);
            auto getXStep(void) const;
            auto getYStep(void) const;
            auto getZStep(void) const;
            auto getXExtent(void) const;
            auto getYExtent(void) const;
            auto getZExtent(void) const;
            /* Field manipulation */
            void setAll(const std::vector<Type>&, const std::vector<Type>&, const std::vector<Type>&);
            void setAll(Type, Type, Type, Type, Type, Type);
            void advectPosition(double, double);

            /* @brief Return the maxmimum x-value of the position field.
               @returns The maximum x-value of the position field.
            */
            Type getXMax()
            {
                return this->xMax;
            }

            /* @brief Return the minimum x-value of the position field.
               @returns The minimum x-value of the position field.
            */
            Type getXMin()
            {
                return this->xMin;
            }

            /* @brief Return the minimum y-value of the position field.
               @returns The minimum y-value of the position field.
            */ 
            Type getYMin()
            {
                return this->yMin;
            }

            /* @brief Return the maximum y-value of the position field.
               @returns The maximum y-value of the position field.
            */
            Type getYMax()
            {
                return this->yMax;
            }

            /* @brief Return the minimum z-value of the position field.
               @returns The minimum z-value of the position field.
            */
            Type getZMin()
            {
                return this->zMin;
            }

            /* @brief Return the maximum z-value of the position field.
               @returns The maximum z-value of the position field.
            */
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

    }; // class
    #include "PositionGetters.hpp"
    #include "Position.inl"
}; // namespace

#endif