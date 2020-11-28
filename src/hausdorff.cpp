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
double hausdorffDistance(std::vector<vectorType>& a, std::vector<vectorType>& b)
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

/*  @brief Computes an array containing whether two strainlines are within some Hausdorff distance H of each other. 
    This is a precursor for determining which strainlines to keep.

    @param[inout] std::vector<std::vector<bool>> relationshipArray Contains whether strainlines are within some tolerance of each other.
    @param[in] 3D array of trajectories - (number of trajectories x steps in trajectory x dimensionality)
    @param[in] tol Tolerance distance for saying two strainlines are 'close'
*/
template <typename vectorType, typename toleranceType>
void computeWithinHausdorffDistance(std::vector<std::vector<bool>>& relationshipArray, std::vector<std::vector<std::vector<vectorType>>>& trajectories, toleranceType& tol)
{
    /* Get size details of the array */
    number_of_trajectories = trajectories.size();
    /* Check the relationship array - is it empty, is it the correct size? */
    if (relationshipArray.empty())
    {
        throw std::runtime_error("The relationship passed to computeWithinHausdorffDistance is empty.");
    } else if (relationshipArray.size() != trajectories.size())
    {
        throw std::runtime_error("The relationship array and trajectories array do not match in size.");
    } else
    {
        for (int relationshipRow = 0; relationshipRow < relationshipArray.size(); ++relationshipRow)
        {
            if (relationshipArray.size() != relationshipArray[relationshipRow].size())
            {
                throw std::runtime_error("The relationship array passed to computeWithinHausdorffDistance is not square.");
            }
        }
    }
    /* OK, now run through every possible trajectory pair in relationship array and see if they are within a given tolerance.
     * First, we assume that all of the values are sufficiently separated, i.e. relationshipArray is all zero/false.
     */
    for (auto relationshipRow : relationshipArray)
    {
        std::fill(relationshipRow.begin(), relationshipRow.end(), false);
    }
    /* Now run through the upper triangular portion of the relationship array (symmetric about diagonal axis) and fill in the 'correct'
    * relationship
    * */
    for (int row = 0; row < relationshipArray.size(); ++row)
    {
        for (int col = row+1; col < relationshipArray.size(); ++col) // Square array
        {
            std::vector<vectorType> trajectoryOnRow = trajectories[row];
            std::vector<vectorType> trajectoryOnCol = trajectories[col];
            double distance = hausdorffDistance(trajectoryOnRow, trajectoryOnCol);
            if (distance < tol)
            {
                relationshipArray[row][col] = true;
            } /* Already initialised relationshipArray to false, so no action to be taken if not true */
        }
    }
    /* Now mirror the operations in the transpose */
    for (int col = 0; col < relationshipArray.size(); ++col)
    {
        for (int row = col + 1; row < relationshipArray.size(); ++row)
        {
            relationshipArray[row][col] = relationshipArray[col][row];
        }
    }
}

#endif