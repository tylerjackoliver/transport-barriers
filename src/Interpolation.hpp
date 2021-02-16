#ifndef __INTERPOLATION_H__
#define __INTERPOLATION_H__

#include <vector>
#include <SPLINTER/definitions.h>
#include <SPLINTER/bsplinebuilder.h>
#include <SPLINTER/bspline.h>
#include <SPLINTER/datatable.h>
#include <iostream>
#include <algorithm>
#include <boost/multi_array.hpp>
#include "Point.hpp"
#include <stdexcept>
#include "Helicity.hpp"
#include "Position.hpp"
#include <chrono>

/* @brief Populate a new boost::multi_array with interpolated helicity values
 * @param[in] HelicityField Boost::Multi_array containing helicity values
 * @param[in] positionField Boost::Multi_array of points containing positions the helicity field is valid at
 * @param[out] reducedHelicityField Boost::Multi_array of points containing positions of the reduced field
 */
template <typename PointType, typename HelicityType>
void createInterpolatedGrid(LCS::Helicity<HelicityType>& helicityField, LCS::position<PointType>& positionField, LCS::position<PointType>& reducedPositionField, LCS::Helicity<HelicityType>& reducedHelicityField)
{
    auto tic = std::chrono::system_clock::now();
    int plane = 0;
    /* Create data table to store samples */
    SPLINTER::DataTable samples;
    /* Iterate through multi_array: add each point/helicity pair to samples */
    int counter = 0;
    for (int row = 0; row < positionField.getXExtent(); ++row)
    {
        for (int col = 0; col < positionField.getYExtent(); ++col)
        {
            SPLINTER::DenseVector sampleToAdd(2);
            Point<PointType> thisPoint = positionField.getValue(row, col, plane);
            HelicityType thisHelicity = helicityField.helicityField_[row][col][plane];
            sampleToAdd(0) = thisPoint.x;
            sampleToAdd(1) = thisPoint.y;
            samples.addSample(sampleToAdd, thisHelicity);
        }
    }
    /* Build the interpolant */
    SPLINTER::BSpline interpolant = SPLINTER::BSpline::Builder(samples).degree(1).build();
    /* Now iterate through the positions in the reduced field and obtain the helicity there */
    for (int row = 0; row < reducedHelicityField.getXExtent(); ++row)
    {
        for (int col = 0; col < reducedHelicityField.getYExtent(); ++col)
        {
            std::vector<PointType> toEvaluate(2);
            Point<PointType> thisPoint = reducedPositionField.getValue(row, col, plane);
            toEvaluate[0] = thisPoint.x;
            toEvaluate[1] = thisPoint.y;
            HelicityType interpolatedHelicity = interpolant.eval(toEvaluate);
            reducedHelicityField.helicityField_[row][col][plane] = interpolatedHelicity;
        }
    }
    auto toc = std::chrono::system_clock::now();
    std::cout << "Completed interpolating the helicity field in " << std::chrono::duration_cast<std::chrono::seconds>(toc - tic).count() << std::endl;
}

#endif
