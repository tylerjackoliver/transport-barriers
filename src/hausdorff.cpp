#ifndef __HAUSDORFF_H__
#define __HAUSDORFF_H__

#include <algorithm>
#include <random>
#include <vector>
#include <stdexcept>
#include <limits>

/* @brief Computes the Hausdorff distance between two vectors a and b which represent trajectories of strainlines.
   @param[inout] a Vector of vector of points
   @param[inout] b Vector of vector of points
   @returns The Hausdorff distance between a and b
*/
template <typename vectorType>
vectorType hausdorffDistance(std::vector<vectorType>& a, std::vector<vectorType>& b)
{
    /* Error checks */
    if (a.empty() || b.empty()) /* Check neither a nor be is empty */
    {
        throw std::runtime_error("One of the vectors passed to hausdorffDistance is empty.");
    } else if ( a[0].size() != b[0].size() ) /* Dimensionality does not match */
    {
        throw std::runtime_error("The vectors passed to hausdorffDistance do not have the same dimensionality.");
    }

    /* Create a random iterator */
    auto seedGenerator = std::random_device {};
    auto rngGenerator = std::default_random_engine { seedGenerator() };
    
    /* Create an array of indices to be shuffled later. Accesing the indices randomly can lead to an earlier
     * break in the distance computation than computing all of the points linearly.
     */
    std::vector<int> indicesA( a.size() );
    std::vector<int> indicesB( b.size() );

    /* Assign values of indices to the vectors created above */
    for (int i = 0; i < indicesA.size(); ++i)
    {
        indicesA[i] = i;
    }
    for (int i = 0; i < indicesB.size(); ++i)
    {
        indicesB[i] = i;
    } /* A and B may not necessarily be of the same length (different number of points in the trajectory), but _will_ have the same dimensionality */

    /* Shuffle the indices arrays */
    std::shuffle(std::begin(indicesA), std::end(indicesB), rngGenerator);
    std::shuffle(std::begin(indicesB), std::end(indicesB), rngGenerator);

    /* Initialise loop variables */
    double cMin, cMax, d;               // Minimum distance, maximum distance, current distance
    bool haveWeBroken = false;          // Have we had a break in the inner loop

    cMin = 0.0;
    cMax = 0.0;

    for (int indexA : indicesA)
    {
        haveWeBroken = false;
        cMin = std::numeric_limits<double>::infinity();  // Set to "infinity" to force update on first loop

        for (int indexB : indicesB)
        {
            d = 0.0; // Reset distance
            /* Get Euclidean distance - keep as squared for now & square at the end to save computation time */
            for (unsigned dimension = 0; dimension < a[0].size(); ++dimension)
            {
                d += (a[indexA][dimension] - b[indexB][dimension]) * (a[indexA][dimension] - b[indexB][dimension]); // Prevent use of std::pow()
            }
            if (d < cMax) // We have an early termination
            {
                haveWeBroken = true;
                break; // Out of inner loop
            }
            /* If we didn't break, set the minimum distance if we've achieved it */
            if (d < cMin) 
            {
                cMin = d;
            }
        }
        if ( isfinite(cMin) && cMin >= cMax && !haveWeBroken ) // We _didn't_ break out of the loop early
        {
            cMax = cMin;
        }
    }
    return std::sqrt(cMax);
}

#endif