#ifndef __HELICITY_H__
#define __HELICITY_H__

#include "Position.hpp"
#include <ostream>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>
#include "Functor.hpp"

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
            void getHelicityBelowThreshold(double, std::vector<std::vector<int>>&);
            void filterSeedPoints(std::vector<std::vector<int>>&);
            void writeSeedPoints(std::vector<std::vector<int>>&);
            void integrateStrainLines(const std::vector<std::vector<int>> &, double);
            
            /* Public data members */
            position <double>& pos_;
            boost::multi_array<double, 3> helicityField_;
            /**/
        private:
            /**/
            void getDominantEigenVector(std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, 
                                        std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, std::vector<Type>&, Eigen::Vector3d&);
            void getHelicity(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, double&);
            void integrate(std::vector<Type>&);
            int sgn(Type&) const;
            void writeStrainline(const std::vector<Point<Type>>&, const int) const;
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

                    double helicity = fabs(curl.dot(eta));
                    this->helicityField_[i][j][k] = helicity;
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
            // std::cout << xLines[i] << std::endl;
            // std::cout << yLines[i] << std::endl;
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
    void Helicity<Type>::integrateStrainLines(const std::vector<std::vector<int>> &indices, double eps)
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

        //#pragma omp parallel for
        for (std::size_t i = 0; i < numStrainlines; ++i)
        {
            /*
            * Define the integrator. Make it a controlled stepper, so the time-step is adjusted if our initial time-step violates
            * initial accuracy constraints.
            */
            typedef std::vector<Type> state_type;
            boost::numeric::odeint::runge_kutta_fehlberg78<state_type> method;
            auto stepper = boost::numeric::odeint::make_controlled(/*reltol*/1e-012, /*abstol*/1e-012, method);
            
            // Extract the initial condition
            std::vector<int> initialIndex = indices[i];
            Point<Type> initialPoint = pos_.getValue(initialIndex[0], initialIndex[1], initialIndex[2]);

            // Convert initial point into state_type
            state_type x0 = {initialPoint.x, initialPoint.y, initialPoint.z};
            state_type previousStep(3);

            // Initialise initial time, current time, time-step 
            double t = 0.0;
            double dt = 0.01; // Guess

            // Instantiate functor class
            functorClass functor;

            // Set problem parameters - initial eigenvector & integration time for eigenvector field
            functor.previousSolution = pos_.getEigenvector(initialIndex[0], initialIndex[1], initialIndex[2]);
            functor.initialTime = pos_.getInitialTime(); // For eigenvector field!
            functor.finalTime = pos_.getFinalTime(); // For eigenvector field!
            functor.xGridSpacing = pos_.getXStep() * pos_.auxiliaryGridSizingFactor;
            functor.yGridSpacing = pos_.getYStep() * pos_.auxiliaryGridSizingFactor;
            functor.zGridSpacing = pos_.getZStep() * pos_.auxiliaryGridSizingFactor;

            // Initialise running helicity
            double helicityTotal = 0.0;
            double helicityAvg = 0.0;
            uint64_t numSteps = 0;   // Step counter - long in case of exceptionally huge strainlines

            // Was the step ok?
            bool success = false;

            // Position history of the particle
            std::vector<Point<double>> strainline;

            // Steps
            double xStep = pos_.getXStep() * pos_.auxiliaryGridSizingFactor;
            double yStep = pos_.getYStep() * pos_.auxiliaryGridSizingFactor;
            double zStep = pos_.getZStep() * pos_.auxiliaryGridSizingFactor;

            // Enter while loop
            while (helicityAvg < eps)
            {
                /* Try and make a step */
                previousStep = x0;
                success = stepper.try_step(functor, x0, t, dt);
                
                // If success is true - IE produced a new x, we need to know the helicity at the new point
                if (success)
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
                    std::vector<Type> leftmZ0(3), upmZ0(3), rightmZ0(3), downmZ0(3), mZmZ0(3), mZ(0);

                    origin = previousStep; origin0 = origin;
                    // Left
                    left = origin; left[0] -= xStep; left0 = left;
                    // Leftleft
                    leftleft = origin; leftleft[0] -= 2. * xStep; leftleft0 = leftleft;
                    // Right
                    right = origin; right[0] += xStep; right0 = right;
                    // Rightright
                    rightright = origin; rightright[0] += 2. * xStep; rightright0 = rightright;
                    // Up
                    up = origin; up[1] += yStep; up0 = up;
                    // Upup
                    upup = origin; up[1] += 2. * yStep; upup0 = upup;
                    // Down
                    down = origin; down[1] -= yStep; down0 = down;
                    // Downdown
                    downdown = origin; downdown[1] -= 2. * yStep; downdown0 = downdown;
                    // Leftup
                    leftup = left; leftup[1] += yStep; leftup0 = leftup;
                    // Leftdown
                    leftdown = left; leftdown[1] -= yStep; leftdown0 = leftdown;
                    // Rightup
                    rightup = right; rightup[1] += yStep; rightup0 = rightup;
                    // Rightdown
                    rightdown = right; rightdown[1] -= yStep; rightdown0 = rightdown;
                    // Left pZ
                    leftpZ = left; leftpZ[2] += zStep; leftpZ0 = leftpZ0;
                    // Left mZ
                    leftmZ = left; leftmZ[2] -= zStep; leftmZ0 = leftmZ;
                    // Right pZ
                    rightpZ = right; rightpZ[2] += zStep; rightpZ0 = rightpZ;
                    // Right mZ
                    rightmZ = right; rightmZ[2] -= zStep; rightmZ0 = rightmZ;
                    // Up pZ
                    uppZ = up; uppZ[2] += zStep; uppZ0 = uppZ;
                    // Up mZ
                    upmZ = up; upmZ[2] -= zStep; upmZ0 = upmZ;
                    // down pZ
                    downpZ = down; downpZ[2] += zStep; downpZ0 = downpZ;
                    // down mZ
                    downmZ = down; downmZ[2] -= zStep; downmZ0 = downmZ;
                    // pZ
                    pZ = origin; pZ[2] += zStep; pZ0 = pZ;
                    // pZpZ
                    pZpZ = pZ; pZ[2] += zStep; pZpZ0 = pZpZ;
                    // mZ
                    mZ = origin; mZ[2] -= zStep; mZ0 = mZ;
                    // mZmZ
                    mZmZ = origin; mZmZ[2] -= zStep; mZmZ0 = mZmZ;

                    // Now we need to integrate each of these points forwards
                    this->integrate(origin); this->integrate(left); this->integrate(leftleft); this->integrate(right); this->integrate(rightright);
                    this->integrate(down); this->integrate(downdown); this->integrate(up); this->integrate(upup);
                    this->integrate(pZ); this->integrate(pZpZ); this->integrate(mZ); this->integrate(mZmZ);

                    this->integrate(leftup); this->integrate(leftdown); this->integrate(rightup); this->integrate(rightdown);
                    this->integrate(leftpZ); this->integrate(rightpZ); this->integrate(downpZ); this->integrate(uppZ);
                    this->integrate(leftmZ); this->integrate(rightmZ); this->integrate(downmZ); this->integrate(upmZ);

                    // Now want to integrate each of the 'x' points forward to get the dominant eigenvector
                    Eigen::Vector3d leftEV, upEV, rightEV, downEV, pZEV, mZEV, originEV;
    
                    // Left
                    this->getDominantEigenVector(&leftleft, &leftup, &origin, &leftdown, &leftpZ, &leftmZ, 
                                                 &leftleft0, &leftup0, &origin0, &leftdown0, &leftpZ0, &leftmZ0, &leftEV);
                    // Up
                    this->getDominantEigenVector(&leftup, &upup, &rightup, &origin, &uppZ, &upmZ, 
                                                 &leftup0, &upup0, &rightup0, &origin0, &uppZ0, &upmZ0, &upEV);
                    // Right
                    this->getDominantEigenVector(&origin, &rightup, &rightright, &rightdown, &rightpZ, &rightmZ, 
                                                 &origin0, &rightup0, &rightright0, &rightdown0, &rightpZ0, &rightmZ0, &rightEV);
                    // Down
                    this->getDominantEigenVector(&leftdown, &origin, &rightdown, &downdown, &downpZ, &downmZ, 
                                                 &leftdown0, &origin0, &rightdown0, &downdown0, &downpZ0, &downmZ0, &downEV);
                    // pZ
                    this->getDominantEigenVector(&leftpZ, &uppZ, &rightpZ, &downpZ, &pZpZ, &origin, 
                                                 &leftpZ0, &uppZ0, &rightpZ0, &downpZ0, &pZpZ0, &origin0, &pZEV);
                    // mZ
                    this->getDominantEigenVector(&leftmZ, &upmZ, &rightmZ, &downmZ, &mZmZ, &origin, 
                                                 &leftmZ0, &upmZ0, &rightmZ0, &downmZ0, &mZmZ0, &origin0, &mZEV);
                    // Origin
                    this->getDominantEigenVector(&left, &up, &right, &down, &pZ, &mZ, &left0, &up0, &right0, &down0, &pZ0, &mZ0, &originEV);

                    // Now get the helicity
                    double helicity;
                    this->getHelicity(&originEV, &leftEV, &upEV, &rightEV, &downEV, &pzEV, &mZEV, &helicity);

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

                    // Update functor class with latest solution
                    functor.previousSolution = x0;

                } // if(success)
            } // while
            this->writeStrainline(&strainline, i);
        } // for
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
    int Helicity<Type>::sgn(Type &in) const
    {
        return (in > 0) - (in < 0);
    }

    template <typename Type>
    void Helicity<Type>::getDominantEigenVector(std::vector<Type>& left, std::vector<Type>& up, std::vector<Type>& right, std::vector<Type>& down, std::vector<Type> &pZ, std::vector<Type> &mZ,
                                std::vector<Type>& left0, std::vector<Type>& up0, std::vector<Type>& right0, std::vector<Type>& down0, std::vector<Type> &pZ0, std::vector<Type> &mZ0, Eigen::Vector3d& ev)
    {
        /* Construct CGST */

        Eigen::Matrix3d<Type> deformation, cgst;

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
    void Helicity<Type>::getHelicity(const Eigen::Vector3d &origin, const Eigen::Vector3d &left, const Eigen::Vector3d &up, const Eigen::Vector3d &right, const Eigen::Vector3d &down, const Eigen::Vector3d &pZ, const Eigen::Vector3d &mZ, double &helicity)
    {
        double xGridSize = pos_.getXStep() * pos_.auxiliaryGridSizingFactor;
        double yGridSize = pos_.getYStep() * pos_.auxiliaryGridSizingFactor;
        double zGridSize = pos_.getZStep() * pos_.auxiliaryGridSizingFactor;

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
    void Helicity<Type>::writeStrainline(const std::vector<Point<Type>>& strainline, const int index) const
    {
        std::ofstream output;
        // Get buffer size required
        std::size_t bufLen = std::snprintf(nullptr, 0, "strainline_%4d", index);
        // Create new character array; using char* allocates on heap in fringe case of huge fnames
        char *buf = new char[bufLen+1];
        std::snprintf(buf, bufLen+1, "strainlines_%4d", i);
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
#endif
