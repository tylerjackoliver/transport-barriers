#ifndef __FTLE_H__
#define __FTLE_H__

#include "Position.hpp"

namespace LCS
{

   template <class Type>
    class FTLE : public field<double>
    {

        private:

            typedef boost::multi_array_types::extent_range extents;

            int sgn(Type val)
            {
                return (val > 0) - (val < 0);
            }

            void integrate(std::vector<Type>& x)
            {
                
                typedef std::vector<Type> state_type;
                int sign = sgn(pos_.getFinalTime() - pos_.getInitialTime());
                boost::numeric::odeint::bulirsch_stoer<state_type> bulirsch_stepper(pos_.getAbsTol(), pos_.getRelTol());
                boost::numeric::odeint::integrate_adaptive(bulirsch_stepper, dynSystem, x, pos_.getInitialTime(), pos_.getFinalTime(), sign*.01, abcFlowObserver);
            }

            void _computeFTLEFiniteDifferencing()
            {

                unsigned prog = 0;
                double coeff = 1.0 / (pos_.getFinalTime() - pos_.getInitialTime());

                #pragma omp parallel for shared(prog)
                for (unsigned i = 0; i < pos_.getXExtent(); ++i)
                {
                    Point<Type> xNext, xPrev, yNext, yPrev, zNext, zPrev;
                    Point<Type> x0Next, x0Prev, y0Next, y0Prev, z0Next, z0Prev;
                    Eigen::Matrix<Type, 3, 3> deformation, cauchy_green;
                    for (unsigned j = 0; j < pos_.getYExtent(); ++j)
                    {
                        for (unsigned k = 0; k < pos_.getZExtent(); ++k){
                            
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
                    std::cout << "Completed processing " << double(prog)/this->nx_ * 100. << "% of the FTLE." << std::endl;
                }
            }

            void _computeFTLEAuxiliaryGrid()
            {
                unsigned prog = 0;
                double coeff = 1.0 / (pos_.getFinalTime() - pos_.getInitialTime());

                #pragma omp parallel for shared(prog)
                for (unsigned i = 0; i < pos_.getXExtent(); ++i)
                {
                    for (unsigned j = 0; j < pos_.getYExtent(); ++j)
                    {
                        for (unsigned k = 0; k < pos_.getZExtent(); ++k){

                            Point<Type> xNext, xPrev, yNext, yPrev, zNext, zPrev;
                            Point<Type> x0Next, x0Prev, y0Next, y0Prev, z0Next, z0Prev;
                            Point<Type> referencePoint, previousReferencePointX, previousReferencePointY, previousReferencePointZ;
                            Eigen::Matrix<Type, 3, 3> deformation, cauchy_green;

                            std::vector<double> reference(3), xPlus(3), xMinus(3), yPlus(3), yMinus(3), zPlus(3), zMinus(3);

                            referencePoint = pos_.getValue(i, j, k);
                            previousReferencePointX = pos_.getValue(i-1, j, k);
                            previousReferencePointY = pos_.getValue(i, j-1, k);
                            previousReferencePointZ = pos_.getValue(i, j, k-1);

                            double xGridSpacing = referencePoint.x - previousReferencePointX.x;
                            double yGridSpacing = referencePoint.y - previousReferencePointY.y;
                            double zGridSpacing = referencePoint.z - previousReferencePointZ.z; 

                            // Get auxiliary grid sizing

                            double auxGridSizeX = xGridSpacing * this->auxiliaryGridSizingFactor;
                            double auxGridSizeY = yGridSpacing * this->auxiliaryGridSizingFactor;
                            double auxGridSizeZ = zGridSpacing * this->auxiliaryGridSizingFactor;

                            reference[0] = referencePoint.x; reference[1] = referencePoint.y;
                            reference[2] = referencePoint.z;

                            xPlus = reference; xPlus[0] += auxGridSizeX; x0Next.vectorToPoint(xPlus);
                            xMinus = reference; xMinus[0] -= auxGridSizeX; x0Prev.vectorToPoint(xMinus);
                            yPlus = reference; yPlus[1] += auxGridSizeY; y0Next.vectorToPoint(yPlus);
                            yMinus = reference; yMinus[1] -= auxGridSizeY; y0Prev.vectorToPoint(yMinus);
                            zPlus = reference; zPlus[2] += auxGridSizeZ; z0Next.vectorToPoint(zPlus);
                            zMinus = reference; zMinus[2] -= auxGridSizeZ; z0Prev.vectorToPoint(zMinus);

                            this->integrate(xPlus); xNext.vectorToPoint(xPlus);  // Implicit casting std::vector -> struct Point
                            this->integrate(xMinus); xPrev.vectorToPoint(xMinus);
                            this->integrate(yPlus); yNext.vectorToPoint(yPlus);
                            this->integrate(yMinus); yPrev.vectorToPoint(yMinus);
                            this->integrate(zPlus); zNext.vectorToPoint(zPlus);
                            this->integrate(zMinus); zPrev.vectorToPoint(zMinus);

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
                    std::cout << "Completed processing " << double(prog)/this->nx_ * 100. << "% of the FTLE." << std::endl;
                }
            }

        public:

            bool useAuxiliaryGrid = 0;
            double auxiliaryGridSizingFactor = 0.25; // % of normal grid size

            FTLE(position<double>& pos): field(pos.getXExtent(), pos.getYExtent(), pos.getZExtent()),
                pos_(pos), ftleField_(boost::extents[range(-1, (int)(pos.getXExtent()+1))][range(-1, (int)(pos.getYExtent()+1))][range(-1, (int)(pos.getZExtent()+1))]){
                };
            
            position<double>& pos_;
            boost::multi_array<double, 3> ftleField_;

            void computeFTLE(){ 
                unsigned xExtent = pos_.getXExtent();
                unsigned yExtent = pos_.getYExtent();
                unsigned zExtent = pos_.getZExtent();
                double coeff = 1./(pos_.getFinalTime() - pos_.getInitialTime());

                std::cout << "Computing the FTLE field...";
                #pragma omp parallel for
                for (unsigned i = 0; i < xExtent; ++i)
                {
                    for (unsigned j = 0; j < yExtent; ++j)
                    {
                        for (unsigned k = 0; k < zExtent; ++k)
                        {
                            // Get the eigenvalue
                            this->ftleField_[i][j][k] = coeff * std::log(pos_.getEigenvalue(i,j,k));
                        }
                    }
                }
                std::cout << "done." << std::endl;
            }

            void writeFTLE()
            {
                std::ofstream output;
                output.open("ftleField.data");

                std::cout << "Writing FTLE values to file...";
                for(unsigned i = 0; i < pos_.getXExtent(); ++i)
                {
                    for(unsigned j = 0; j < pos_.getYExtent(); ++j)
                    {
                        for(unsigned k = 0; k < pos_.getZExtent(); ++k)
                        {
                            output << ftleField_[i][j][k] << std::endl;
                        }
                    }
                }
                std::cout << "done.";
                output.close();
            }
    }; // Class
}; // Namespace

#endif