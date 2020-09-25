#include "domain.hpp"
#include <boost/numeric/odeint.hpp>

domain::domain(unsigned int _xDim, unsigned int _yDim, unsigned int _zDim) : xDim(_xDim), yDim(_yDim), zDim(_zDim)
{
    initialDomain.xDim = xDim; 
    initialDomain.yDim = yDim; 
    initialDomain.zDim = zDim;
    initialDomain.reallocateMemory();

    finalDomain.xDim = xDim;
    finalDomain.yDim = yDim;
    finalDomain.zDim = zDim;
    finalDomain.reallocateMemory();
};

// We want to ensure maximum throughput. Thus, construct a few different version depending on
// whether OpenMPI is used or not

#if defined(_OPENMP)

void domain::initialisePosition(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
    #pragma omp parallel for shared(xmax, xmin)
    for (int i = 0; i <= xDim+2; ++i)
    {
        this->initialDomain.x[i] = (xmax-xmin) * (i-1) / this->initialDomain.xDim;
    }

    #pragma omp parallel for
    for (int i = 0; i < yDim + 2; ++i)
    {
        this->initialDomain.y[i] = (ymax-ymin) * (i-1) / this->initialDomain.yDim;
    }

    #pragma omp parallel for
    for (int i = 0; i < zDim + 2; ++i)
    {
        this->initialDomain.z[i] = (zmax-zmin) * (i-1) / this->initialDomain.zDim;
    }

}

#else

// Non-OpenMP versions
void domain::initialisePosition(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
    #pragma simd
    for (unsigned i = 0; i < this->xDim+2; ++i)
    {
        this->initialDomain.x[i] = (xmax-xmin) * (i-1) / this->xDim;
    }

    #pragma simd
    for (unsigned i = 0; i < this->yDim+2 + 2; ++i)
    {
        this->initialDomain.y[i] = (ymax-ymin) * (i-1) / this->yDim;
    }

    #pragma simd
    for (unsigned i = 0; i < this->zDim+2 + 2; ++i)
    {
        this->initialDomain.z[i] = (zmax-zmin) * (i-1) / this->zDim;
    }
}

#endif

// template <doublename double>
void domain::initialisePosition(std::vector<double> range)
{
    double xmin = range[0];
    double xmax = range[1];
    double ymin = range[2];
    double ymax = range[3];
    double zmin = range[4];
    double zmax = range[5];

    this->initialisePosition(xmin, xmax, ymin, ymax, zmin, zmax);
}

domain::~domain()
{};

domain::advect(double absTol, double relTol, double t0, double tf)
{
    #pragma omp parallel for
    for (int i = 0; i < (this->xDim + 2); ++i)
    {
        for (int j = 0; j < (this->yDim + 2); ++j)
        {
            for (int k = 0; k < (this->zDim + 2); ++k)
            {
                double stateVector[dimensionality];
                // Reversing the order of the class storage may be more effective - we shall see
                x[0] = this->initialDomain.x[i];
                x[1] = this->initialDomain.y[j];
                x[2] = this->initialDomain.z[k];

                // Get the sign of the integration
                int sign = ( (tf-t0) > 0 ) - ( (tf-t0) < 0);

                // Advect forward
                boost::numeric::odeint::bulirsch_stoer<(double *)> bulirsch_stepper(absTol, relTol);
                boost::numeric::odeint::integrate_adaptive(bulirsch_stepper, dynSystem, stateVector, t0, tf, sign);

                this->finalDomain.x[i]


            }
        }
    }
}