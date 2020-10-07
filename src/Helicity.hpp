#ifndef __HELICITY_H__
#define __HELICITY_H__

#include "Position.hpp"
#include <ostream>

namespace LCS
{

    template <class Type>
    class Helicity
    {

        public:
            /* @brief Class constructor for the Helicity class.
            @param[in] pos Reference to a populated position class
            */
            Helicity(position<double>& _pos) : pos_(_pos)
            {
                // Correctly initialise the boost::multiarray as public class member
                typedef boost::multi_array<double, 3> tempType;
                typename tempType::extent_gen extents;
                this->helicityField_.resize(extents[pos_.getXExtent()+2][pos_.getYExtent()+2][pos_.getZExtent()+2]);
                this->helicityField_.reindex(-1);
            }
            /* Public method signatures */

            void computeHelicityField();
            void writeHelicityField() const;
            void writeHelicityField(std::string);
            void readHelicityField();
            void readHelicityField(std::string);
            void writeHelicityBelowThreshold(double) const;
            
            /* Public data members */

            position <double>& pos_;
            boost::multi_array<double, 3> helicityField_;
    };
    template <typename Type>
    void Helicity<Type>::computeHelicityField()
    {
        std::cout << "Computing helicity field...";

        #pragma omp parallel for
        for (unsigned i = 0; i < pos_.getXExtent(); ++i)
        {
            for (unsigned j = 0; j < pos_.getYExtent(); ++j)
            {
                for (unsigned k = 0; k < pos_.getZExtent(); ++k)
                {
                    Eigen::Vector3d eta(3), curl(3), xNext, xPrev, yNext, yPrev, zNext, zPrev;
                    Point<Type> x0Next, x0Prev, y0Next, y0Prev, z0Next, z0Prev;

                    eta = pos_.getEigenvector(i,j,k);

                    xNext = pos_.getEigenvector(i+1, j, k);
                    xPrev = pos_.getEigenvector(i-1, j, k);

                    yNext = pos_.getEigenvector(i, j+1, k);
                    yPrev = pos_.getEigenvector(i, j-1, k);

                    zNext = pos_.getEigenvector(i, j, k+1);
                    zPrev = pos_.getEigenvector(i, j, k-1);

                    x0Next = pos_.getValue(i+1, j, k);
                    x0Prev = pos_.getValue(i-1, j, k);

                    y0Next = pos_.getValue(i, j+1, k);
                    y0Prev = pos_.getValue(i, j-1, k);

                    z0Next = pos_.getValue(i, j, k+1);
                    z0Prev = pos_.getValue(i, j, k-1);

                    double dfzdy = (yNext(2) - yPrev(2)) / (y0Next.y - y0Prev.y);
                    double dfydz = (zNext(1) - zPrev(1)) / (z0Next.z - z0Prev.z);

                    double dfxdz = (zNext(0) - zPrev(0)) / (z0Next.z - z0Prev.z);
                    double dfzdx = (xNext(2) - xPrev(2)) / (x0Next.x - x0Prev.x);

                    double dfydx = (xNext(1) - xPrev(1)) / (x0Next.x - x0Prev.x);
                    double dfxdy = (yNext(0) - yPrev(0)) / (y0Next.y - y0Prev.y);

                    curl(0) = (dfzdy - dfydz);
                    curl(1) = (dfxdz - dfzdx);
                    curl(2) = (dfydx - dfxdy);

                    this->helicityField_[i][j][k] = fabs(curl.dot(eta));
                }
            }
        }

        std::cout << "done." << std::endl;
    }

    template <typename Type>
    void Helicity<Type>::writeHelicityField() const
    {
        std::ofstream output;
        output.open("../helicity.data");
        for (unsigned i = 0; i < pos_.getXExtent(); ++i)
        {
            for (unsigned j = 0; j < pos_.getYExtent(); ++j)
            {
                for (unsigned k = 0; k < pos_.getZExtent(); ++k)
                {
                    output << this->helicityField_[i][j][k] << std::endl;
                }
            }
        }

        output.close();
    }

    template <typename Type>
    void Helicity<Type>::writeHelicityBelowThreshold(double eps) const
    {
        std::ofstream output;
        output.open("../helicityMask.data");
        for (unsigned i = 0; i < pos_.getXExtent(); ++i)
        {
            for (unsigned j = 0; j < pos_.getYExtent(); ++j)
            {
                for (unsigned k = 0; k < pos_.getZExtent(); ++k)
                {
                    double helicity = this->helicityField_[i][j][k];
                    if (helicity < eps) output << 1 << std::endl;
                    else output << 0 << std::endl;
                }
            }
        }
        output.close();
    }
};
#endif
