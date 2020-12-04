#ifndef __INTERPOLATION_H__
#define __INTERPOLATION_H__

#include <vector>
#include <SPLINTER/bspline.h>
#include <SPLINTER/datatable.h>
#include <iostream>
#include <algorithm>
#include <boost/multi_array.hpp>
#include "Point.hpp"
#include <stdexcept>
/* @brief Populate a new boost::multi_array with interpolated helicity values
 * @param[in] HelicityField Boost::Multi_array containing helicity values
 * @param[in] positionField Boost::Multi_array of points containing positions the helicity field is valid at
 * @param[out] reducedHelicityField Boost::Multi_array of points containing positions of the reduced field
 */
template <typename PointType, typename HelicityType>
void createInterpolatedGrid(boost::multi_array<HelicityType, 2>& helicityField, boost::multi_array<Point<PointType>, 2>& positionField, boost::multi_array<HelicityType, 2>& reducedHelicityField)
{
    /* Create data table to store samples */
    SPLINTER::DataTable samples;
    /* Ensure helicityField and positionField are of the same size */
    if ( helicityField.shape() != positionField.shape() )
    {
        throw std::runtime_error("The helicity field and position field are not of the same size in createInterpolator");
    }
    /* Iterate through multi_array: add each point/helicity pair to samples */
    for (int row = 0; row < positionField.shape()[0]; ++row)
    {
        for (int col = 0; col < positionField.shape()[1]; ++col)
        {
            SPLINTER::DenseVector sampleToAdd(3);
            Point<PointType> thisPoint = positionField[row][col];
            HelicityType thisHelicity = helicityField[row][col];
            sampleToAdd(0) = thisPoint.x;
            sampleToAdd(1) = thisPoint.y;
            sampleToAdd(2) = thisPoint.z;
            samples.addSample(sampleToAdd, thisHelicity);
        }
    }
    /* Build the interpolant */
    SPLINTER::BSpline interpolant = SPLINTER::BSpline::Builder(samples).degree(3).build();
    
    /* Now iterate through the positions in the reduced field and obtain the helicity there */
    for (int row = 0; row < reducedHelicityField.shape()[0]; ++row)
    {
        for (int col = 0; col < reducedHelicityField.shape()[1]; ++col)
        {
            std::vector<PointType> toEvaluate(3);
            Point<PointType> thisPoint = reducedHelicityField[row][col];
            toEvaluate[0] = thisPoint.x;
            toEvaluate[1] = thisPoint.y;
            toEvaluate[2] = thisPoint.z;
            HelicityType interpolatedHelicity = interpolant.eval(toEvaluate);
            reducedHelicityField[row][col] = interpolatedHelicity;
        }
    }
}

#endif
