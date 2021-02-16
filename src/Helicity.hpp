#ifndef __HELICITY_H__
#define __HELICITY_H__

#include "Position.hpp"
#include <ostream>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>
#include "Functor.hpp"
#include <unistd.h>
#include <chrono>
#include "Frechet.hpp"

namespace LCS
{
    /* @brief Contains a series of methods for computing quantities relating to the Helicity field for the flow.
    Helicity is computed with respect to a reference position field, and can be computed either from pre-determined eigenvalues, or
    integrated manually as part of the class.*/
    template <class Type>
    class Helicity
    {

        public:
            /* @brief Class constructor for the Helicity class.
            @param[in] pos Reference to a populated position class
            */
            Helicity(position<double>& _pos) : pos_(_pos)
            {
                /* Initialise boost::multi_array with the correct extents and rebase to -1 for boundary conditions */
                typedef boost::multi_array<double, 3> tempType;
                typename tempType::extent_gen extents;
                this->helicityField_.resize(extents[pos_.getXExtent()+2][pos_.getYExtent()+2][pos_.getZExtent()+2]);
                this->helicityField_.reindex(-1);
                /* Initialise the integration steps for computing auxiliary grids based off of the position field */
                _xStepIntegration = pos_.getXStep() * pos_.auxiliaryGridSizingFactor;
                _yStepIntegration = pos_.getYStep() * pos_.auxiliaryGridSizingFactor;
                _zStepIntegration = pos_.getZStep() * pos_.auxiliaryGridSizingFactor;
            }
            /* Public method signatures */
            void computeHelicityField();
            void writeHelicityField() const;
            void writeHelicityField(std::string);
            void readHelicityField();
            void computeHelicityFieldAuxiliaryGrid();
            void readHelicityField(std::string);
            void writeHelicityBelowThreshold(double) const;
            void getHelicityBelowThreshold(double, std::vector<std::vector<int>>&);
            void filterSeedPoints(std::vector<std::vector<int>>&);
            void writeSeedPoints(std::vector<std::vector<int>>&);
            void integrateStrainLines(const std::vector<std::vector<int>>&, double, std::vector<std::vector<Point<double>>>&);
            int getXExtent(void);
            int getYExtent(void);
            int getZExtent(void);
            void setXStepIntegration(double);
            void setYStepIntegration(double);
            void setZStepIntegration(double);
            /* Public data members */
            position <double>& pos_;
            boost::multi_array<double, 3> helicityField_;
            /**/
        private:
            /* Private member signatures */
            void getDominantEigenVector(std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, 
                                        std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, Eigen::Vector3d&);
            void getDominantEigenVector(std::vector<Type>&, Eigen::Vector3d&);
            void getHelicity(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, double&);
            void integrate(std::vector<Type>&);
            int sgn(Type) const;
            void writeStrainline(const std::vector<Point<Type>>&, const int, const int) const;
            double _xStepIntegration;
            double _yStepIntegration;
            double _zStepIntegration;
    };

    /* @brief Setter function for the Helicity integration. 
       @param[in] xStep The step to use in the auxiliary grid.
    */
    template <typename Type>
    void Helicity<Type>::setXStepIntegration(double xStep)
    {
        this->_xStepIntegration = xStep;
    }

    /* @brief Setter function for the y-step to be used in constructing the auxiliary grid.
       @param[in] yStep The step to use in the y-direction for the auxiliary grid
    */
    template <typename Type>
    void Helicity<Type>::setYStepIntegration(double yStep)
    {
        this->_yStepIntegration = yStep;
    }

    /* @brief Setter function for the z-step to be used in constructing the auxiliary grid.
       @param[in] zStep The step to use in the z-direction for the auxiliary grid.
    */
    template <typename Type>
    void Helicity<Type>::setZStepIntegration(double zStep)
    {
        this->_zStepIntegration = zStep;
    }

    /* @brief Getter function for the x-extent of the domain of the helicity field.
       @returns The x-extent of the domain field.*/
    template <typename Type>
    int Helicity<Type>::getXExtent()
    {
        return pos_.getXExtent();
    }

    /* @brief Getter function for the y-extent of the domain of the helicity field.
       @returns The y-extent of the domain field.*/
    template <typename Type>
    int Helicity<Type>::getYExtent()
    {
        return pos_.getYExtent();   
    }

    /* @brief Getter function for the z-extent of the domain of the helicity field.
       @returns The z-extent of the domain field.*/
    template <typename Type>
    int Helicity<Type>::getZExtent()
    {
        return pos_.getZExtent();
    }

    /* @brief Computes the helicity field over the domain by computing an auxiliary grid of eigenvectors about each point. Thus, there are actually
     (6 x 6) integrations of the trajectory performed per point, since each auxiliary grid requires _it's own_ auxiliary grid. Not quite as intensive as mining bitcoin,
     but not far off.
    */
    template <typename Type>
    void Helicity<Type>::computeHelicityFieldAuxiliaryGrid()
    {
        std::cout << "Computing helicity field...";
        std::cout.flush();

        int numCompleted = 0;
        #pragma omp parallel for shared(numCompleted)
        for (unsigned i = 0; i < pos_.getXExtent(); ++i)
        {
            for (unsigned j = 0; j < pos_.getYExtent(); ++j)
            {
                for (unsigned k = 0; k < pos_.getZExtent(); ++k)
                {
                    Eigen::Vector3d eta(3), curl(3), xNext, xPrev, yNext, yPrev, zNext, zPrev;
                    Point<Type> x0, x0Next, x0Prev, y0Next, y0Prev, z0Next, z0Prev;

                    /* Construct original grid */
                    x0 = pos_.getValue(i, j, k);
                    x0Next = x0; x0Next.x += this->_xStepIntegration;
                    x0Prev = x0; x0Prev.x -= this->_xStepIntegration;
                    y0Next = x0; y0Next.y += this->_yStepIntegration;
                    y0Prev = x0; y0Prev.y -= this->_yStepIntegration;
                    z0Next = x0; z0Next.z += this->_zStepIntegration;
                    z0Prev = x0; z0Prev.z -= this->_zStepIntegration; 

                    /* Compute the eigenvector about each of these points. This also uses an auxiliary grid, so does result in unnecessary computation of overlapping trajectories,
                       but computational effort is not currently an issue. */
                    Helpers::makeGridGetEigenVector(x0, this->_xStepIntegration, this->_yStepIntegration, this->_zStepIntegration, eta, 
                                                    this->pos_.getAbsTol(), this->pos_.getRelTol(), this->pos_.getInitialTime(), this->pos_.getFinalTime());
                    Helpers::makeGridGetEigenVector(x0Next, this->_xStepIntegration, this->_yStepIntegration, this->_zStepIntegration, xNext, 
                                                    this->pos_.getAbsTol(), this->pos_.getRelTol(), this->pos_.getInitialTime(), this->pos_.getFinalTime());
                    Helpers::makeGridGetEigenVector(x0Prev, this->_xStepIntegration, this->_yStepIntegration, this->_zStepIntegration, xPrev, 
                                                    this->pos_.getAbsTol(), this->pos_.getRelTol(), this->pos_.getInitialTime(), this->pos_.getFinalTime());
                    Helpers::makeGridGetEigenVector(y0Next, this->_xStepIntegration, this->_yStepIntegration, this->_zStepIntegration, yNext, 
                                                    this->pos_.getAbsTol(), this->pos_.getRelTol(), this->pos_.getInitialTime(), this->pos_.getFinalTime());
                    Helpers::makeGridGetEigenVector(y0Prev, this->_xStepIntegration, this->_yStepIntegration, this->_zStepIntegration, yPrev, 
                                                    this->pos_.getAbsTol(), this->pos_.getRelTol(), this->pos_.getInitialTime(), this->pos_.getFinalTime());
                    Helpers::makeGridGetEigenVector(z0Next, this->_xStepIntegration, this->_yStepIntegration, this->_zStepIntegration, zNext, 
                                                    this->pos_.getAbsTol(), this->pos_.getRelTol(), this->pos_.getInitialTime(), this->pos_.getFinalTime());
                    Helpers::makeGridGetEigenVector(z0Prev, this->_xStepIntegration, this->_yStepIntegration, this->_zStepIntegration, zPrev, 
                                                    this->pos_.getAbsTol(), this->pos_.getRelTol(), this->pos_.getInitialTime(), this->pos_.getFinalTime());

                    /* Compute curl */
                    double dfzdy = (yNext(2) - yPrev(2)) / (y0Next.y - y0Prev.y);
                    double dfydz = (zNext(1) - zPrev(1)) / (z0Next.z - z0Prev.z);

                    double dfxdz = (zNext(0) - zPrev(0)) / (z0Next.z - z0Prev.z);
                    double dfzdx = (xNext(2) - xPrev(2)) / (x0Next.x - x0Prev.x);

                    double dfydx = (xNext(1) - xPrev(1)) / (x0Next.x - x0Prev.x);
                    double dfxdy = (yNext(0) - yPrev(0)) / (y0Next.y - y0Prev.y);

                    curl(0) = (dfzdy - dfydz);
                    curl(1) = (dfxdz - dfzdx);
                    curl(2) = (dfydx - dfxdy);

                    double helicity = fabs(curl.dot(eta));
                    this->helicityField_[i][j][k] = helicity;
                    pos_.eigenvectorField_[i][j][k] = eta;
                }
            }
            #pragma omp atomic 
            numCompleted++;
            std::cout << "Completed " << static_cast<double>(numCompleted) / pos_.getXExtent()* 100 << "% of the helicity." << std::endl;
        }
        std::cout << "done." << std::endl;
    }

    /* @brief Compute the helicity field using the eigenvectors already computed in the position field.
    */
    template <typename Type>
    void Helicity<Type>::computeHelicityField()
    {
        std::cout << "Computing helicity field...";

        /* We only compute the interior points of the domain here due to boundary conditions, so can go from zero -> extent. */
        #pragma omp parallel for
        for (unsigned i = 0; i < pos_.getXExtent(); ++i)
        {
            for (unsigned j = 0; j < pos_.getYExtent(); ++j)
            {
                for (unsigned k = 0; k < pos_.getZExtent(); ++k)
                {
                    Eigen::Vector3d eta(3), curl(3), xNext, xPrev, yNext, yPrev, zNext, zPrev;
                    Point<Type> x0Next, x0Prev, y0Next, y0Prev, z0Next, z0Prev;

                    /* Construct the eigenvector grid about the central point */
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

                    /* Compute curl, as above. */
                    double dfzdy = (yNext(2) - yPrev(2)) / (y0Next.y - y0Prev.y);
                    double dfydz = (zNext(1) - zPrev(1)) / (z0Next.z - z0Prev.z);

                    double dfxdz = (zNext(0) - zPrev(0)) / (z0Next.z - z0Prev.z);
                    double dfzdx = (xNext(2) - xPrev(2)) / (x0Next.x - x0Prev.x);

                    double dfydx = (xNext(1) - xPrev(1)) / (x0Next.x - x0Prev.x);
                    double dfxdy = (yNext(0) - yPrev(0)) / (y0Next.y - y0Prev.y);

                    curl(0) = (dfzdy - dfydz);
                    curl(1) = (dfxdz - dfzdx);
                    curl(2) = (dfydx - dfxdy);

                    double helicity = fabs(curl.dot(eta));
                    this->helicityField_[i][j][k] = helicity;
                }
            }
        }
        std::cout << "done." << std::endl;
    }

    /* @brief Write the computed helicity field to an output file. by default, this is 'helicity.data' */
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

    template <typename Type>
    void Helicity<Type>::getHelicityBelowThreshold(double eps, std::vector<std::vector<int>> &points)
    {
        for (int i = 0; i < pos_.getXExtent(); ++i)
        {
            for (int j = 0; j < pos_.getYExtent(); ++j)
            {
                for (int k = 0; k < pos_.getZExtent(); ++k)
                {
                    if (this->helicityField_[i][j][k] < eps) 
                    {
                        std::vector<int> indices{i, j, k};
                        points.push_back(indices);
                    }
                }
            }
        }
    }

    template <typename Type>
    void Helicity<Type>::filterSeedPoints(std::vector<std::vector<int>> &points)
    {
        std::vector<std::vector<int>> newPoints;

        // Get x/y coordinates of lines to check
        Point<Type> min = pos_.getValue(0, 0, 0);
        Point<Type> max = pos_.getValue(pos_.getXExtent()-1, pos_.getYExtent()-1, pos_.getZExtent()-1);

        Type minX = min.x; Type minY = min.y;
        Type maxX = max.x; Type maxY = max.y;

        std::vector<Type> xLines(4);
        std::vector<Type> yLines(4);

        for (int i = 0; i < xLines.size(); ++i)
        {
            xLines[i] = (i+1) * (maxX-minX) / (xLines.size()+1);
            yLines[i] = (i+1) * (maxY-minY) / (yLines.size()+1);
        }

        // Now check intersections with the lines
        for (int i = 0; i < points.size(); ++i)
        {
            int ii = points[i][0];
            int jj = points[i][1];
            int kk = points[i][2];

            // Check if point intersects either equi-spaced x- or y-lines
            Point<Type> samplePoint = this->pos_.getValue(ii, jj, kk);
            Type x = samplePoint.x;
            Type y = samplePoint.y;
            Type z = samplePoint.z;

            // Check intersection with x
            [&]{ // LABEL 1
                for (int j = 0; j < xLines.size(); ++j)
                {
                    if (fabs(x - xLines[j]) < 1e-02)
                    {
                        newPoints.push_back(points[i]);
                        return; // Out of Lambda (LABEL 1)
                    } else if (fabs(y-yLines[j]) < 1e-02)
                    {
                        newPoints.push_back(points[i]);
                        return; // Out of Lambda (LABEL 1)
                    }
                }
            }();
        }
        points = newPoints;
    }

    template <typename Type>
    void Helicity<Type>::writeSeedPoints(std::vector<std::vector<int>>& indices)
    {
        std::ofstream output;
        output.open("seedPoints.data");
        // Now check intersections with the lines
        for (int i = 0; i < indices.size(); ++i)
        {
            int ii = indices[i][0];
            int jj = indices[i][1];
            int kk = indices[i][2];

            // Check if point intersects either equi-spaced x- or y-lines
            Point<Type> samplePoint = this->pos_.getValue(ii, jj, kk);
            Type x = samplePoint.x;
            Type y = samplePoint.y;
            Type z = samplePoint.z;

            output << x << std::endl;
            output << y << std::endl;
            output << z << std::endl;
        }
    }

    template <typename Type>
    void Helicity<Type>::integrateStrainLines(const std::vector<std::vector<int>> &indices, double eps, std::vector<std::vector<Point<double>>>& allStrainlines)
    {
        /*
         * We're going to go through every point in indices and integrate the strainlines corresponding to the
         * initial seed point.
         * 
         * We need a numerical integrator (with step size control) to integrate the tensorline very carefully. A
         * custom class is used with overloaded operator() to also memorise the past solution of the ODE, as is required.
         * 
         * The integration will stop when:
         *  - The running average of helicity (taken w.r.t. the number of function evaluations) goes over the tolerance epsilon
         */

        /* Enter into the loop */
        std::size_t numStrainlines = indices.size();
        unsigned int strainLinesCompleted = 0;

        std::cout << "~~~~~~~~~~~~~~~~~~~ STRAINLINE INTEGRATION ~~~~~~~~~~~~~~~~~" << std::endl;
        std::cout << "There are " << numStrainlines << " strainlines to integrate." << std::endl;

        auto timeStart = std::chrono::high_resolution_clock::now();
        // std::vector<std::vector<Point<double>>> allStrainlines;
        #pragma omp parallel for shared(strainLinesCompleted)
        for (std::size_t i = 0; i < numStrainlines; ++i)
        {
            std::vector<std::vector<Point<double>>> bothDirections;

            // Extract the initial condition
            std::vector<int> initialIndex = indices[i];
            Point<Type> initialPoint = pos_.getValue(initialIndex[0], initialIndex[1], initialIndex[2]);

            /* Has this point come close to an already-computed strainline?
                Need while rather than for in case size has increased between iterations
            */
            bool isOK = false;
            int numStrainlinesTested = 0;
            while (numStrainlinesTested < allStrainlines.size())
            {
                // isOK = Helpers::distanceBelowThreshold( allStrainlines[numStrainlinesTested], initialPoint, 10. );
                // if (isOK) break;
                numStrainlinesTested++;
            }
            isOK = false;
            if (isOK && strainLinesCompleted != 0)
            {
                std::cout << "skipped " << strainLinesCompleted << std::endl;
                #pragma omp critical (strainlineupdate)
                {
                    strainLinesCompleted++;
                }
                continue;
            }
            for (int direction = -1; direction < 2; direction+=2)
            {
                /*
                * Define the integrator. Make it a controlled stepper, so the time-step is adjusted if our initial time-step violates
                * initial accuracy constraints.
                */
                typedef std::vector<double> state_type;
                boost::numeric::odeint::runge_kutta_dopri5<state_type> method;
                auto stepper = boost::numeric::odeint::make_controlled(/*reltol*/1e-013, /*abstol*/1e-013, method);
                boost::numeric::odeint::runge_kutta_dopri5<state_type> step2;

                // Convert initial point into state_type
                state_type x0 = {initialPoint.x, initialPoint.y, initialPoint.z};
                state_type previousStep(3);

                // Initialise initial time, current time, time-step 
                double t = 0.0;
                double dt = 0.001; // Guess

                // Instantiate functor class
                functorClass functor;

                // Set problem parameters - initial eigenvector & integration time for eigenvector field
                Eigen::Vector3d initNormal, initUnitNormal;
                Point<Type> initpZ = pos_.getValue(initialIndex[0], initialIndex[1], initialIndex[2] + 1);
                Point<Type> initmZ = pos_.getValue(initialIndex[0], initialIndex[1], initialIndex[2] - 1);
            
                initNormal[0] = initpZ.x - initmZ.x;
                initNormal[1] = initpZ.y - initmZ.y;
                initNormal[2] = initpZ.z - initmZ.z;

                initUnitNormal = initNormal.normalized();

                functor.previousSolution = direction * initUnitNormal.cross(pos_.getEigenvector(initialIndex[0], initialIndex[1], initialIndex[2]));
                functor.initialTime = pos_.getInitialTime(); // For eigenvector field!
                functor.finalTime = pos_.getFinalTime(); // For eigenvector field!

                // Steps
                double xStep = this->_xStepIntegration;
                double yStep = this->_yStepIntegration;
                double zStep = this->_zStepIntegration;

                functor.xGridSpacing = xStep;
                functor.yGridSpacing = yStep;
                functor.zGridSpacing = zStep;

                // Initialise running helicity
                double helicityTotal = this->helicityField_[initialIndex[0]][initialIndex[1]][initialIndex[2]];
                double helicityAvg = helicityTotal;
                double pi = 4.0 * std::atan(1.0);
                double length = 0.0;
                unsigned long numSteps = 1;
                const unsigned long long maxSteps = 250e6;

                // Was the step ok?
                boost::numeric::odeint::controlled_step_result result;

                // Position history of the particle
                std::vector<Point<double>> strainline;

                // Enter while loop
                while (helicityAvg < eps)
                {
                    /* Try and make a step */
                    previousStep = x0;
                    result = stepper.try_step(functor, x0, t, dt);
                    dt = std::min(dt, 0.01);
                
                    // If success is true - IE produced a new x, we need to know the helicity at the new point
                    if (result == boost::numeric::odeint::success)
                    {
                        /*
                        * We now need to compute the gradient of the helicity field. This requires eigenvectors of the surrounding 4 points, and we therefore
                        * need to perform a total of 13 integrations. If an eigenvector point is x, and a 'normal' integral point is '.', the structure is:
                        *                              up
                        *                              .
                        *                          .   x   .
                        *           left       .   x   o   x   .    right
                        *                          .   x   .
                        *                              .
                        *                            down
                        * With o the origin point. The central four points also need to be duplicated in the +Z and -Z axes to get the 3D terms.
                        * 
                        * We denote each of the following vectors by their location in this 'diamond' - e.g. leftleft being the leftmost, leftup being the top-left point, etc.
                        * There is also a pZ (+Z) and mZ (-Z) direction.
                        * 
                        * I am so sorry for what you're about to read. - JT 10.10.2020
                        */
                        std::vector<Type> origin(3), left(3), leftleft(3), up(3), upup(3), down(3), downdown(3), right(3), rightright(3);
                        std::vector<Type> origin0(3), left0(3), leftleft0(3), up0(3), upup0(3), down0(3), downdown0(3), right0(3), rightright0(3);

                        std::vector<Type> leftup(3), leftdown(3), rightup(3), rightdown(3);
                        std::vector<Type> leftup0(3), leftdown0(3), rightup0(3), rightdown0(3);

                        std::vector<Type> leftpZ(3), uppZ(3), rightpZ(3), downpZ(3), pZpZ(3), pZ(3);
                        std::vector<Type> leftpZ0(3), uppZ0(3), rightpZ0(3), downpZ0(3), pZpZ0(3), pZ0(3);

                        std::vector<Type> leftmZ(3), upmZ(3), rightmZ(3), downmZ(3), mZmZ(3), mZ(3);
                        std::vector<Type> leftmZ0(3), upmZ0(3), rightmZ0(3), downmZ0(3), mZmZ0(3), mZ0(3);

                        // Now want to integrate each of the 'x' points forward to get the dominant eigenvector
                        Eigen::Vector3d leftEV, upEV, rightEV, downEV, pZEV, mZEV, originEV;

                        // Up
                        origin = previousStep;
                        up = origin; up[1] += yStep; this->getDominantEigenVector(up, upEV);
                        down = origin; down[1] -= yStep; this->getDominantEigenVector(down, downEV);
                        left = origin; left[0] -= xStep; this->getDominantEigenVector(left, leftEV);
                        right = origin; right[0] += xStep; this->getDominantEigenVector(right, rightEV);
                        pZ = origin; pZ[2] += zStep; 
                        mZ = origin; mZ[2] -= zStep;

                        // Update functor class with latest solution+
                        Eigen::Vector3d normal, unitNormal;
                        for (unsigned j = 0; j < pZ.size(); ++j)
                        {
                            normal[j] = pZ[j] - mZ[j];
                        }
                        unitNormal = normal.normalized();
                        this->getDominantEigenVector(mZ, mZEV);
                        this->getDominantEigenVector(pZ, pZEV);
                        this->getDominantEigenVector(origin, originEV);

                        // Now get the helicity
                        double helicity;
                        this->getHelicity(originEV, leftEV, upEV, rightEV, downEV, pZEV, mZEV, helicity);
                        // Update the running average
                        numSteps++;
                        helicityTotal += helicity;
                        helicityAvg = helicityTotal / static_cast<double>(numSteps);

                        // Add the updated position to our history
                        Point<Type> latestPoint;
                        latestPoint.x = x0[0];
                        latestPoint.y = x0[1];
                        latestPoint.z = x0[2];
                        strainline.push_back(latestPoint);

                        Eigen::Vector3d newSol = unitNormal.cross(originEV);
                        double innerProduct = functor.previousSolution.dot(newSol);
                        std::vector<double> diffVec(3);
                        std::transform(x0.begin(), x0.end(), previousStep.begin(), diffVec.begin(), std::minus<double>());
                        length += std::sqrt(diffVec[0] * diffVec[0] + diffVec[1]*diffVec[1] + diffVec[2]*diffVec[2]);
                        functor.previousSolution = sgn(innerProduct) * newSol; // Flip direction if required

                        // Fix maximum number of steps
                        if (numSteps > maxSteps) helicityAvg = 100000000.;
                    } 
                }
                bothDirections.push_back(strainline);   // Add to the pair
                std::cout << "Completed integrating strainline " << strainLinesCompleted << " of " << numStrainlines << " after " << numSteps << "steps." << std::endl;
            } //direction 
            /* Reverse the direction of the first trajectory to get one smooth continuous line (and remove repeated point) */
            std::reverse(bothDirections[1].begin(), bothDirections[1].end());
            bothDirections[1].pop_back();
            /* Create new vector containing the entire trajectory */
            std::vector<Point<double>> entireTrajectory;
            entireTrajectory.insert(entireTrajectory.begin(), bothDirections[1].begin(), bothDirections[1].end());
            entireTrajectory.insert(entireTrajectory.end(), bothDirections[0].begin(), bothDirections[0].end());
            #pragma omp critical (strainlineupdate)
            {
                strainLinesCompleted++;
                this->writeStrainline(entireTrajectory, strainLinesCompleted, 1);
                allStrainlines.push_back(entireTrajectory);
            }
        } // for
        auto timeEnd = std::chrono::high_resolution_clock::now();
        std::cout << "Integrating " << numStrainlines << " strainlines complete. Time required: " << std::chrono::duration_cast<std::chrono::minutes>(timeEnd-timeStart).count() << " minutes." << std::endl;
    } // class member function

    template <typename Type>
    void Helicity<Type>::integrate(std::vector<Type> &x)
    {
        typedef std::vector<Type> state_type;
        int sign = this->sgn(pos_.getFinalTime() - pos_.getInitialTime());
        boost::numeric::odeint::bulirsch_stoer<state_type> bulirsch_stepper(pos_.getAbsTol(), pos_.getRelTol());
        boost::numeric::odeint::integrate_adaptive(bulirsch_stepper, dynSystem, x, pos_.getInitialTime(), pos_.getFinalTime(), sign*.01, abcFlowObserver);
    }

    template <typename Type>
    int Helicity<Type>::sgn(Type in) const
    {
        return (in > 0) - (in < 0);
    }

    template <typename Type>
    void Helicity<Type>::getDominantEigenVector(std::vector<Type>& left, std::vector<Type>& up, std::vector<Type>& right, std::vector<Type>& down, std::vector<Type> &pZ, std::vector<Type> &mZ,
                                std::vector<Type>& left0, std::vector<Type>& up0, std::vector<Type>& right0, std::vector<Type>& down0, std::vector<Type> &pZ0, std::vector<Type> &mZ0, Eigen::Vector3d& ev)
    {
        /* Construct CGST */

        Eigen::Matrix3d deformation, cauchy_green;

        // Deformation tensor
        deformation(0,0) = (right[0]-left[0]) / (right0[0]-left0[0]);
        deformation(0,1) = (up[0]-down[0]) / (up0[1]-down0[1]);
        deformation(0,2) = (pZ[0]-mZ[0]) / (pZ0[2]-mZ0[2]);
        
        deformation(1,0) = (right[1]-left[1]) / (right0[0]-left0[0]);
        deformation(1,1) = (up[1]-down[1]) / (up0[1]-down0[1]);
        deformation(1,2) = (pZ[1]-mZ[1]) / (pZ0[2]-mZ0[2]);

        deformation(2,0) = (right[2]-left[2]) / (right0[0]-left0[0]);
        deformation(2,1) = (up[2]-down[2]) / (up0[1]-down0[1]);
        deformation(2,2) = (pZ[2]-mZ[2]) / (pZ0[2]-mZ0[2]);
                            
        cauchy_green = deformation.transpose() * deformation;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(cauchy_green);
        ev = solver.eigenvectors().col(2).real();
    }

    template <typename Type>
    void Helicity<Type>::getDominantEigenVector(std::vector<Type>& x0, Eigen::Vector3d& ev)
    {
        Eigen::Matrix<Type, 3, 3> deformation, cauchy_green;
        std::vector<double> xPlus(3), xMinus(3), yPlus(3), yMinus(3), zPlus(3), zMinus(3);

        double auxGridSizeX = pos_.getXStep() * pos_.auxiliaryGridSizingFactor;
        double auxGridSizeY = pos_.getYStep() * pos_.auxiliaryGridSizingFactor;
        double auxGridSizeZ = pos_.getZStep() * pos_.auxiliaryGridSizingFactor;

        xPlus = x0; xPlus[0] += auxGridSizeX;
        xMinus = x0; xMinus[0] -= auxGridSizeX;
        yPlus = x0; yPlus[1] += auxGridSizeY;
        yMinus = x0; yMinus[1] -= auxGridSizeY;
        zPlus = x0; zPlus[2] += auxGridSizeZ;
        zMinus = x0; zMinus[2] -= auxGridSizeZ;

        this->integrate(xPlus);
        this->integrate(xMinus);
        this->integrate(yPlus);
        this->integrate(yMinus);
        this->integrate(zPlus);
        this->integrate(zMinus);

        // deformation tensor 
        deformation(0,0) = (xPlus[0]-xMinus[0]) / (2. * auxGridSizeX);
        deformation(0,1) = (yPlus[0]-yMinus[0]) / (2. * auxGridSizeY);
        deformation(0,2) = (zPlus[0]-zMinus[0]) / (2. * auxGridSizeZ);
        
        deformation(1,0) = (xPlus[1]-xMinus[1]) / (2. * auxGridSizeX);
        deformation(1,1) = (yPlus[1]-yMinus[1]) / (2. * auxGridSizeY);
        deformation(1,2) = (zPlus[1]-zMinus[1]) / (2. * auxGridSizeZ);

        deformation(2,0) = (xPlus[2]-xMinus[2]) / (2. * auxGridSizeX);
        deformation(2,1) = (yPlus[2]-yMinus[2]) / (2. * auxGridSizeY);
        deformation(2,2) = (zPlus[2]-zMinus[2]) / (2. * auxGridSizeZ);
        
        cauchy_green = deformation.transpose() * deformation;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(cauchy_green);
        ev = solver.eigenvectors().col(2).real();
    }

    template <typename Type>
    void Helicity<Type>::getHelicity(const Eigen::Vector3d &origin, const Eigen::Vector3d &left, const Eigen::Vector3d &up, const Eigen::Vector3d &right, const Eigen::Vector3d &down, const Eigen::Vector3d &pZ, const Eigen::Vector3d &mZ, double &helicity)
    {
        Eigen::Vector3d curl;

        double xGridSize = pos_.getXStep();
        double yGridSize = pos_.getYStep();
        double zGridSize = pos_.getZStep();

        double dfzdy = (up(2) - down(2)) / (2. * yGridSize);
        double dfydz = (pZ(1) - mZ(1)) / (2. * zGridSize);

        double dfxdz = (pZ(0) - mZ(0)) / (2. * zGridSize);
        double dfzdx = (right(2) - left(2)) / (2. * xGridSize);

        double dfydx = (right(1) - left(1)) / (2. * xGridSize);
        double dfxdy = (up(0) - down(0)) / (2. * yGridSize);

        curl(0) = (dfzdy - dfydz);
        curl(1) = (dfxdz - dfzdx);
        curl(2) = (dfydx - dfxdy);

        helicity = fabs(curl.dot(origin));
    }

    template <typename Type>
    void Helicity<Type>::writeStrainline(const std::vector<Point<Type>>& strainline, const int index, const int direction) const
    {
        std::ofstream output;
        // Get buffer size required
        std::size_t bufLen = std::snprintf(nullptr, 0, "../strainlines/strainlines_%d", index);
        // Create new character array; using char* allocates on heap in fringe case of huge fnames
        char *buf = new char[bufLen+1];
        snprintf(buf, bufLen+1, "../strainlines/strainlines_%d", index);
        // Open file
        output.open(buf);
        // Get vector size
        std::size_t vectorLen = strainline.size();
        // Write to output
        for (std::size_t i = 0; i < vectorLen; ++i)
        {
            Point <Type> tmp = strainline[i];
            output << tmp.x << std::endl;
            output << tmp.y << std::endl;
            output << tmp.z << std::endl;
        }
        // Clean
        output.close();
        delete[] buf;
    }

}; // namespace

template <typename Type>
std::ostream& operator<<(std::vector<Type>& in, const std::ostream& out)
{
    out << "(";
    for (unsigned i = 0; i < in.size(); ++i) out << in[i] << ",";
    return out << ")";
}

#endif
