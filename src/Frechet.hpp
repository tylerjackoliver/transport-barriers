#ifndef __FRECHET_H__
#define __FRECHET_H__

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <chrono>
#include "Point.hpp"

namespace Frechet
{
    template <typename T>
    T distanceMetric(const Point<T>& a, const Point<T>& b)
    {
        T sum = 0.0;
        sum += (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
        return std::sqrt(sum);
    }

    /* Computes the optimized distance matrix as in Devogele, T., Esnault, M., Etienne, L., & Lardy, F. (2017) */
    template <typename T>
    T computeDistanceMatrix( std::vector<Point<T>>& l1, std::vector<Point<T>>& l2, std::vector<std::vector<T>>& distanceMatrix)
    {
        /* Counters and sentinels */
        int i = 0, j = 0;
        T diagMax = 0.0;
        /* Compute trajectory lengths, remainders and quotients */
        int n = l1.size();
        int m = l2.size();
        if (n < m) /* Need to swap the trajectories */
        {
            std::swap(n, m);
            std::swap(l1, l2);
        }
        int q = static_cast<int>(n / m);
        int r = n % m;
        /* Resize the distance matrix to be of the right size. */
        distanceMatrix.reserve(n);
        for (int idx = 0; idx < n; ++idx) distanceMatrix.push_back(std::vector<T>(m));
        /* Now construct the 'almost diagonal' */
        for (i = 0; i <= r * (q + 1); ++i)
        {
            j = i / (q + 1);
            distanceMatrix[i][j] = distanceMetric(l1[i], l2[j]);
            diagMax = std::max(distanceMatrix[i][j], diagMax);
        }
        /* Construct the end of the diagonal */
        for (i = r * (q + 1); i <= (n-1); ++i)
        {
            j = (i-r) / q;
            distanceMatrix[i][j] = distanceMetric(l1[i], l2[j]);
            diagMax = std::max(distanceMatrix[i][j], diagMax);
        }
        /* Construct the upper-right matrix */
        int jMin = 0;
        for (i = 0; i <= (n-1); ++i)
        {
            j = i;
            T d;
            while ( (j < m || d <= diagMax) && j < jMin)
            {
                d = distanceMetric(l1[i], l2[j]);
                if (d < diagMax)
                {
                    distanceMatrix[i][j] = d;
                    j++;
                }
            };
            jMin = j;
        }
        /* Lower left matrix */
        int iMin = 0;
        for (j = 0; j <= (m-1); ++j)
        {
            i = j;
            T d = 0.0;
            do 
            {
                d = distanceMetric(l1[i], l2[j]);
                if (d < diagMax)
                {
                    distanceMatrix[i][j] = d;
                    i++;
                }
            } while ( (i < n || d <= diagMax) && (i < iMin) ) ;
            iMin = i;
        }
        return diagMax;
    }


    /* @brief Compute the Frechet matrix as in Devogele, T., Esnault, M., Etienne, L., & Lardy, F. (2017). Optimized Discrete Fréchet Distance between trajectories. */
    template <typename T>
    void computeFrechetMatrix(std::vector<Point<T>>& l1, std::vector<Point<T>>& l2, std::vector<std::vector<T>>& frechetMatrix)
    {
        int i = 0, j = 0, jMin = 0;
        /* Compute the distance matrix - it will swap l1 and l2 if required,
        * so get n and m after this in case there's been a swap. */
        std::vector<std::vector<T>> distanceMatrix;
        T diagMax = computeDistanceMatrix(l1, l2, distanceMatrix);
        int n = l1.size(), m = l2.size();
        /* Resize the Frechet Matrix */
        frechetMatrix.reserve(n);
        for (unsigned idx = 0; idx < n; ++idx) frechetMatrix.push_back(std::vector<T>(m));
        frechetMatrix[0][0] = distanceMatrix[0][0];
        for (i = 1; i <= (n-1); ++i)
        {
            j = jMin;
            while ( distanceMatrix[i][j] == 0 ) j++;
            jMin = j;
            while ( distanceMatrix[i][j] != 0 && j < m)
            {
                // T minimum = numeric_limits<T>::max() - 1; // Set to 'infinity'
                T minimum = 1000.0;
                if (i > 0 && j > 0 && distanceMatrix[i-1][j-1] != 0)
                {
                    minimum = frechetMatrix[i-1][j-1];
                }
                if (i > 0 && distanceMatrix[i-1][j] != 0)
                {
                    minimum = std::min(minimum, frechetMatrix[i-1][j]);
                }
                if (j > 0 && distanceMatrix[i][j-1] != 0)
                {
                    minimum = std::min(minimum, frechetMatrix[i][j-1]);
                }
                frechetMatrix[i][j] = std::max(minimum, distanceMatrix[i][j]);
                j++;
            }
        }
    }

    /* @brief Computes the Frechet distance between two trajectories. */
    template <typename T>
    T frechetDistance(std::vector<Point<T>>& l1, std::vector<Point<T>>& l2)
    {
        std::vector<std::vector<T>> frechetMatrix;
        computeFrechetMatrix(l1, l2, frechetMatrix);
        return frechetMatrix[ l1.size() - 1 ][ l2.size() - 1];
    }

    /* @brief Filter the strainlines found by Frechet distance.
     * @param[in] std::vector<std::vector<Point<T>>> containing all the strainlines.
     * @param[in] minDistance minimum distance that two trajectories must be separated by to be 'independent'
     * @param[out] filteredStrainlines std::vector<std::vector<Point<T>>> containing the filtered strainlines.
     */
    template <typename T>
    void filterOnFrechet(std::vector<std::vector<Point<T>>>& allStrainlines, double minDistance, std::vector<std::vector<Point<T>>>& filteredStrainlines)
    {
        /* OK, so have a list of all the possible strainlines in the dataset. We're going to go through one-by-one and check whether each trajectory
           is similar to another one we've already seen. If it is, we arbitarily remove the other trajectory and keep the one we were examining.
           We'll need to initialise the array with the first strainline, so make sure it's a good one!
         
            - Would taking the longest strainline be best? Should in theory result in the same answer, but something to try
        */
        /* Check there are actually strainlines in teh dataset. */
        if ( !allStrainlines.size() ) throw std::runtime_error("Sorry, there were no strainlines to filter!");
        filteredStrainlines.push_back(allStrainlines[0]);
        std::cout << "Got here " << std::endl;
        /* There is a slight chance that we get similar trajectories included twice, depending on racing. However, that's a price to pay here for the parallelism
        speedup. */
        // #pragma omp parallel for shared(allStrainlines, filteredStrainlines) schedule(dynamic, 4) // Want worksharing to be able to reassign when number of trajectories added gets low. 
        for (int traj = 1; traj < ( allStrainlines.size() - 1 ); ++traj) // -1 as last trajectory will have already been checked when we get there
        {
            std::vector<Point<T>> thisTrajectory = allStrainlines[traj];
            /* By default, assumine it's going into the new strainline database. */
            bool toInclude = true;
            for (int toCompare = traj + 1; toCompare < allStrainlines.size(); ++toCompare)
            {
                double distance = frechetDistance(thisTrajectory, allStrainlines[toCompare]);
                /* If this trajectory is within another that's already in the database, stop computing, don't add this one. */
                if (distance < minDistance)
                {
                    toInclude = false;
                    break;
                }
            }
            if (toInclude) filteredStrainlines.push_back(thisTrajectory);
            std::cout << "Completed comparing trajectory " << traj << " of " << allStrainlines.size() << '\n';
        }
    };

    /* Writes all the filtered strainlines to files */
    template <typename T>
    void writeFrechet(std::vector<std::vector<Point<T>>>& allStrainlines)
    {
        for (int i = 0; i < allStrainlines.size(); ++i)
        {
            std::ofstream output;
            output.open(std::string("../filtered_strainlines/strainline_") + std::to_string(i));
            for (auto& step : allStrainlines[i])
            {
                output << step.x << ",";
                output << step.y << ",";
                output << step.z << ",";
                output << "\n";
            }
            output.close();
            std::cout << "Completed writing strainline " << i << " of " << allStrainlines.size() << '\n';
        }
    }
};

#endif
