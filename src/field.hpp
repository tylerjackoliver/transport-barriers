#ifndef __FIELD_H__
#define __FIELD_H__

#include <omp.h>
#include <boost/multi_array.hpp>
#include <eigen3/Eigen/Eigenvalues>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include "forceFunction.hpp"
#include <fstream>
#include "Point.hpp"
#include <iostream>

/* template <typename Type>
void advectPositionGPUDriver(field_<Type>& initPos, field_<Type>& updatedPos, double initTime, double finalTime, \
double absTol, double relTol);
*/

template <typename Type>
void operator<<(std::ostream& out, std::vector<Type>& val)
{
    std::cout << "(";
    for(unsigned i=0; i<val.size(); ++i){std::cout << val[i] << ", ";}
    std::cout << ")";
}

namespace LCS
{
    template <class Type>
    class field
    {

        private:
                    
        protected:

            typedef std::vector<Type> state_type;           // Integrator state type
            unsigned nx_;                                   // X-extent
            unsigned ny_;                                   // Y-extent
            unsigned nz_;                                   // Z-extent
            typedef boost::multi_array<Point<Type>, 3> tempType;
            tempType data_;       // Three-dimensional array of points
            typedef boost::multi_array_types::extent_range range;

        public:

            /* @brief Constructor for the field class.
             *
             * @param nx The number of points in x
             * @param ny The number of points in y
             * @param nz The number of points in z
             */
            field(unsigned nx, unsigned ny, unsigned nz) : nx_(nx), ny_(ny), nz_(nz)
            {
                // Initialise data to be the correct size based on xDim/yDim/zDim

                typename tempType::extent_gen extents;
                this->data_.resize(extents[nx+2][ny+2][nz+2]);

                // Set data to be -1-indexed
                this->data_.reindex(-1);
            }

            /* @brief Retrieve elements of the domain
             *
             * @param i The index of the point in x
             * @param j The index of the point in y
             * @param k The index of the point in z
             */
            Point<Type> getValue(unsigned i, unsigned j, unsigned k)
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
            void setValue(Point<Type> value, unsigned i, unsigned j, unsigned k)
            {

                data_[i][j][k] = value;

            }

    };

}; // namespace

#endif