#ifndef __POSITION_H__
#define __POSITION_H__

#include "field.hpp"

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

            template <typename T> int sgn(T val)
            {
                return (val > 0) - (val < 0);              
            }

        public:

            boost::multi_array<Point<Type>, 3> updatedPosition_;
            double initTime_;
            double finalTime_;
            double absTol_ = 1e-08; // Integrator settings
            double relTol_ = 1e-08; // Integrator settings

            /* Constructor for initialising the position field
            * @param[in] nx The number of points in the \f$x\f$-direction
            * @param[in] ny The number of points in the \f$y\f$-direction
            * @param[in] nz The number of points in the \f$z\f$-direction
            */

            position(unsigned nx, unsigned ny, unsigned nz): field<Type>(nx, ny, nz){
            
                typedef boost::multi_array<Point<Type>, 3> tempType;
                typename tempType::extent_gen extents;
                this->updatedPosition_.resize(extents[nx+2][ny+2][nz+2]);
                updatedPosition_.reindex(-1);

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
                
                std::generate(xrange.begin(), xrange.end(), \
                        [&]{ return xmin + (i++) * (xmax-xmin) / (this->nx_-1);});
                std::generate(yrange.begin(), yrange.end(), \
                        [&]{ return ymin + (j++) * (ymax-ymin) / (this->ny_-1);});

                // Correction for if we're on a single z-plane (planar flow)

                if (this->nz_ < 2)
                {
                    Type gap = (zmax - zmin) / (this->nz_+1);
                    std::generate(zrange.begin(), zrange.end(), 
                        [&]{ return zmin + (k++) * gap;});
                } else
                {
                    std::generate(zrange.begin(), zrange.end(), \
                        [&]{ return zmin + (k++) * (zmax-zmin) / (this->nz_-1);});
                }

                setAll(xrange, yrange, zrange);
                    
            }

            /* @brief Set the initial time associated with the particles of the field
            *
            * @param[in] initTime Initial time of the field
            */
            void setInitialTime(double initTime){

                this->initTime_ = initTime;

            }

            /* @brief Set the final time associated with the advection of the particles of the field.
            *
            * @param[in] finalTime Final time associated with the advection of the particles of the field.
            */
            void setFinalTime(double finalTime){

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
                    std::cout << "Completed " << double(prog)/(this->nx_+1) * 100. << "%." << std::endl;

                } // i

            } // void advectPosition()

            void advectPositionGPU(double absTol, double relTol){

                advectPositionGPUDriver(this->data_, updatedPosition_, initTime_, finalTime_, absTol, relTol);
                Type temp = updatedPosition_[3][3][3];
            
            }

            /* @brief Get the number of points in the 'x'-direction
            *
            * @returns The number of points in the 'x'-direction
            */
            auto getXExtent()
            {

                return this->nx_;

            }

            /* @brief Get the number of points in the 'y'-direction
            *
            * @returns The number of points in the 'y'-direction
            */
            auto getYExtent()
            {

                return this->ny_;

            }

            /* @brief Get the number of points in the 'z'-direction
            *
            * @returns The number of points in the 'z'-direction
            */
            auto getZExtent()
            {

                return this->nz_;

            }

            /* @brief Write the updated position to a file "field.data"
            */
            void writeToFile()
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
            void writeToFile(std::string fName)
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
                std::ofstream input;
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
                std::ofstream input;
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

        protected:


    }; // class

}; // namespace
 

#endif
