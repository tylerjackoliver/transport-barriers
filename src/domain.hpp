#ifndef __DOMAIN_H__
#define __DOMAIN_H__

#include "Point_array.hpp"
#include <iostream>

// template <class double>
class domain
{
    public:

        domain(unsigned int xDim, unsigned int yDim, unsigned int zDim);
        ~domain();

        void initialisePosition(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
        void initialisePosition(std::vector<double> range);
    
        Points<double> initialDomain;
        Points<double> finalDomain;

    private:

        unsigned int xDim;
        unsigned int yDim;
        unsigned int zDim;

};

#include "domain.tpp"

#endif