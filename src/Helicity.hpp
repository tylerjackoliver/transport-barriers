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
            void getDominantEigenVector(std::vector<Type>&, Eigen::Vector3d&);
            void getHelicity(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, double&);
            void integrate(std::vector<Type>&);
            int sgn(Type) const;
            void writeStrainline(const std::vector<Point<Type, 3>>&, const int, const int) const;
    };

    template <typename Type>
    void Helicity<Type>::computeHelicityField()
    {
        std::cout << "Computing helicity field...";
        unsigned long progress = 0;
        #pragma omp parallel for shared(progress)
        for (unsigned i = 0; i < pos_.getXExtent(); ++i)
        {
            for (unsigned j = 0; j < pos_.getYExtent(); ++j)
            {
                for (unsigned k = 0; k < pos_.getZExtent(); ++k)
                {
                    /* Obtain auxiliary grid step sizes. */
                    double xStep = pos_.getXStep() * pos_.auxiliaryGridSizingFactor;
                    double yStep = pos_.getYStep() * pos_.auxiliaryGridSizingFactor;
                    double zStep = pos_.getZStep() * pos_.auxiliaryGridSizingFactor;

                    /* Construct auxiliary grid */
                    Eigen::Vector3d thisEV, xNext, xPrev, yNext, yPrev, zNext, zPrev;
                    Point<Type, 3> x0Next, x0Prev, y0Next, y0Prev, z0Next, z0Prev, thisPoint;

                    thisPoint = pos_.getValue(i, j, k);
                    x0Next = thisPoint; x0Next.state[0] += xStep;
                    x0Prev = thisPoint; x0Prev.state[0] -= xStep;
                    y0Next = thisPoint; y0Next.state[1] += yStep;
                    y0Prev = thisPoint; y0Prev.state[1] -= yStep;
                    z0Next = thisPoint; z0Next.state[2] += zStep;
                    z0Prev = thisPoint; z0Prev.state[2] -= zStep;

                    /* Obtain eigenvector at each of these points */
                    this->getDominantEigenVector(thisPoint, thisEV)
                    this->getDominantEigenVector(x0Next, xNext);
                    this->getDominantEigenVector(x0Prev, xPrev);
                    this->getDominantEigenVector(y0Next, yNext);
                    this->getDominantEigenVector(y0Prev, yPrev);
                    this->getDominantEigenVector(z0Next, zNext);
                    this->getDominantEigenVector(z0Prev, zPrev);

                    /* Now get the helicity using this Point. */
                    double helicity;
                    this->getHelicity(thisEV, xPrev, yNext, xNext, yPrev, zNext, zPrev, helicity);
                    this->helicityField_[i][j][k] = helicity;
                }
            }
            progress++;
            std::cout << "Completed computing " << dynamic_cast<double>(progress) / pos_.getXExtent() * 100 << "% of the helicity field." << '\n';
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
        Point<Type, 3> min = pos_.getValue(0, 0, 0);
        Point<Type, 3> max = pos_.getValue(pos_.getXExtent()-1, pos_.getYExtent()-1, pos_.getZExtent()-1);

        Type minX = min.state[0]; Type minY = min.state[1];
        Type maxX = max.state[0]; Type maxY = max.state[1];

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
            Point<Type, 3> samplePoint = this->pos_.getValue(ii, jj, kk);
            Type x = samplePoint.state[0];
            Type y = samplePoint.state[1];
            Type z = samplePoint.state[2];

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
            Point<Type, 3> samplePoint = this->pos_.getValue(ii, jj, kk);
            Type x = samplePoint.state[0];
            Type y = samplePoint.state[1];
            Type z = samplePoint.state[2];

            output << x << std::endl;
            output << y << std::endl;
            output << z << std::endl;
        }
    }

    /* @brief Integrate the strainlines given in indices for as long as the running average of Helicity remains below a tolerance epsilon.
       @param[in] relTol Template parameter; defaults to 1e-013. Relative tolerance of the numerical integration.
       @param[in] absTol Template parameter; defaults to 1e-013. Absolute tolerance of the numerical integration.
       @param[in] eps Threshold to use for the helicity computation.
    */
    template <typename Type, double relTol = 1e-013, double absTol=1e-013>
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
        unsigned int strainLinesCompleted = 0;

        std::cout << "~~~~~~~~~~~~~~~~~~~ STRAINLINE INTEGRATION ~~~~~~~~~~~~~~~~~" << std::endl;
        std::cout << "There are " << numStrainlines << " strainlines to integrate." << std::endl;

        auto timeStart = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for shared(strainLinesCompleted)
        for (std::size_t i = 0; i < numStrainlines; ++i)
        {
            /* Store the result of the integration in both directions. */
            std::vector<std::vector<Point<double>>> bothDirections;
            /* Start by going negative, then go poisitive. */
            for (int direction = -1; direction < 2; direction+=2)
            {
                /*
                * Define the integrator. Make it a controlled stepper, so the time-step is adjusted if our initial time-step violates
                * initial accuracy constraints.
                */
                typedef std::vector<double> state_type;
                boost::numeric::odeint::runge_kutta_dopri5<state_type> method;
                auto stepper = boost::numeric::odeint::make_controlled(relTol, absTol, method);
                boost::numeric::odeint::runge_kutta_dopri5<state_type> step2;
                
                /* Extract the initial condition from the indices array */
                std::vector<int> initialIndex = indices[i];
                state_type x0 = pos_.getValue(initialIndex[0], initialIndex[1], initialIndex[2]).state;

                /* Convert the initial point from a Point object t*/
                state_type x0 = initialPoint.state;
                state_type previousStep(3);

                /* Initialise initial time and initial guess for the time-step. */
                double t = 0.0;
                double dt = 0.001;

                /* Instantiate a functor class, which defines the RHS of the ODE for the integrator to use. It also provides a method for storing
                 * the current unit normal to the plane, which is a requirement of the force function for projection onto the correct plane.
                 */
                functorClass functor;

                /* Set theunit normal in functor. */
                Eigen::Vector3d initNormal, initUnitNormal;
                /* +- Z */
                Point<Type, 3> initpZ = pos_.getValue(initialIndex[0], initialIndex[1], initialIndex[2] + 1);
                Point<Type, 3> initmZ = pos_.getValue(initialIndex[0], initialIndex[1], initialIndex[2] - 1);
            
                /* Construct dimensional normal. */
                initNormal[0] = initpZ.state[0] - initmZ.state[0];
                initNormal[1] = initpZ.state[1] - initmZ.state[1];
                initNormal[2] = initpZ.state[2] - initmZ.state[2];

                /* Construct unit normal. */
                initUnitNormal = initNormal.normalized();

                /* Initialise the value of the 'previous solution' in Functor. This is used to keep the strainline integration globally smooth - to prevent reversals,
                 * the value of Functor at each step will be compared to the value at the next step.
                 */
                functor.previousSolution = direction * initUnitNormal.cross(pos_.getEigenvector(initialIndex[0], initialIndex[1], initialIndex[2]));
                /* Set the initial and final time for the nested eigenvector determination integration in the Functor */
                functor.initialTime = pos_.getInitialTime(); 
                functor.finalTime = pos_.getFinalTime(); 

                /* Set step sizes for use in the Functor nested integration and in the helicity determination later. */
                double xStep = pos_.getXStep() * pos_.auxiliaryGridSizingFactor;
                double yStep = pos_.getYStep() * pos_.auxiliaryGridSizingFactor;
                double zStep = pos_.getZStep() * pos_.auxiliaryGridSizingFactor;

                functor.xGridSpacing = xStep;
                functor.yGridSpacing = yStep;
                functor.zGridSpacing = zStep;

                /* Initialise total and running helicity. */
                double helicityTotal = this->helicityField_[initialIndex[0]][initialIndex[1]][initialIndex[2]];
                double helicityAvg = helicityTotal;
                /* Define pi for use later in terminating the integration at a maximum strainline length. */
                double pi = 4.0 * std::atan(1.0);
                /* Track the length of the transfer as part of one of the stopping conditions. */
                double length = 0.0;
                /* Track the number of integration steps taken as part of the stopping conditions. */
                unsigned long numSteps = 1;
                /* Set the maximum number of steps the integrator can take in computing the strainline. */
                const unsigned long long maxSteps = 250e6;

                /* Status variable for checking if steps completed successfully. */
                boost::numeric::odeint::controlled_step_result result;

                /* Position history of the particle. */
                std::vector<Point<double, 3>> strainline;

                /* Integrate as long as the running average of the helicity is below the threshold. */
                while (helicityAvg < eps)
                {
                    /* Try and make a step */
                    previousStep = x0;
                    result = stepper.try_step(functor, x0, t, dt);
                    /* Constrain the step size to prevent the integrator 'skipping over' important features. */
                    dt = std::min(dt, 0.01);
                
                    /* If we were properly able to make a step, we need to compute the new value of helicity at the given time. */
                    if (result == boost::numeric::odeint::success)
                    {
                        /* 
                         * Now we need to determine the helicity at this new point. Well do this in a similar (read: exactly the same) way to the computeHelicity() function:
                         * we'll create a new auxiliary grid of points and call getDominantEigenVector() at each of these points to compute the helicity. 
                         * 
                         * Nomenclature : x0 is the solution at the new point we want to test.
                         * 
                         * Note: this could be made more efficient as we're essentially recomputing trajectories at the overlapping grid points, but leave as-is for now.
                         * 
                         */
                        state_type xNext = x0; xNext[0] += xStep; // Initialise to reference point, add on the grid sizing to translate to the correct position.
                        state_type xPrev = x0; xPrev[0] -= xStep;
                        state_type yNext = x0; yNext[1] += yStep;
                        state_type yPrev = x0; xPrev[1] -= yStep;
                        state_type zNext = x0; zNext[2] += zStep;
                        state_type zPrev = x0; zPrev[2] -= zStep;

                        /* Get the dominant eigenvectors at each of these points. */
                        Eigen::Matrix<Type, 3, 1> x0EV, xNextEV, xPrevEV, yNextEV, yPrevEV, zNextEV, zPrevEV;

                        this->getDominantEigenVector(x0, x0EV);
                        this->getDominantEigenVector(xNext, xNextEV);
                        this->getDominantEigenVector(xPrev, xPrevEV);
                        this->getDominantEigenVector(yNext, yNextEV);
                        this->getDominantEigenVector(yPrev, yPrevEV);
                        this->getDominantEigenVector(zNext, zNextEV);
                        this->getDominantEigenVector(zPrev, zPrevEV);

                        /* Now compute the helicity */
                        double helicity;
                        this->getHelicity(x0, xPrev, yNext, xNext, yPrev, zNext, zPrev, helicity);
                        /* Update the number of steps, get the running average of helicity */
                        numSteps++;
                        helicityTotal += helicity;
                        helicityAvg = helicityTotal / dynamic_cast<double>(numSteps); // Prevent integer division (extreme but best to make sure!)

                        /* Doing quick (superfluous) checks on the helicity average now can save us having to do much more extra work later on if we're over the threshold with this latest step.*/
                        if (helicityAvg > eps) continue;
                        if (numSteps > maxSteps)
                        {
                            helicityAvg = 1.01 * eps; // Warning if eps is set ridiculously huge..!
                        }
                        
                        /* Get the new unit normal at this point */
                        Eigen::Matrix<Type, 3, 1> unitNormal;
                        for (unsigned j = 0; j < zNext.size(); ++j)
                        {
                           unitNormal[j] = zNext[j] - zPrev[j];
                        }
                        unitNormal = unitNormal.normalized();

                        /* Add the updated position to our history. */
                        Point<Type, 3> latestPoint(x0);
                        strainline.push_back(latestPoint);

                        /* Update previousSolution in the Functor with this new step */
                        Eigen::Matrix<double, 3, 1> newSolution = unitNormal.cross(x0EV);
                        double innerProduct = functor.previousSolution.dot(newSolution);
                        functor.previousSolution = sgn(innerProduct) * newSolution;

                        /* Compute the new length of the strainline. */
                        std::vector<double> differenceVector(3);
                        std::transform( x0.begin(), x0.end(), previousStep.begin(), diffVec.begin(), std::minus<double>() );
                        length += std:;sqrt(diffVec[0] * diffVec[0] + diffVec[1] * diffVec[1] + diffVec[2] * diffVec[2]);

                        /* Check if the length is greater than the maximum length */
                        if (length >= 10 * pi) helicityAvg = 100000000.;
                    }
                }
                /* Add this direction of the strainline to the total strainline trajectory. */
                bothDirections.push_back(strainline);
            }
            /* Reverse the direction of the first trajectory to get one smooth continuous line (and remove repeated point) */
            std::reverse(bothDirections[1].begin(), bothDirections[1].end());
            bothDirections[1].pop_back();

            /* Create new vector containing the entire trajectory */
            std::vector<Point<double>> entireTrajectory;
            entireTrajectory.insert(entireTrajectory.begin(), bothDirections[1].begin(), bothDirections[1].end());
            entireTrajectory.insert(entireTrajectory.end(), bothDirections[0].begin(), bothDirections[0].end());

            #pragma omp critical
            {
                strainLinesCompleted++;
                this->writeStrainline(entireTrajectory, strainLinesCompleted, 1);
            }

            std::cout << "Completed integrating strainline " << strainLinesCompleted << " of " << numStrainlines << std::endl;
        }
        auto timeEnd = std::chrono::high_resolution_clock::now();
        std::cout << "Integrating " << numStrainlines << " strainlines complete. Time required: " << std::chrono::duration_cast<std::chrono::minutes>(timeEnd-timeStart).count() << " minutes." << std::endl;
    }

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
    void Helicity<Type>::writeStrainline(const std::vector<Point<Type, 3>>& strainline, const int index, const int direction) const
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
            output << tmp.state[0] << std::endl;
            output << tmp.state[1] << std::endl;
            output << tmp.state[2] << std::endl;
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
