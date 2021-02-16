#include <iostream>
#define EIGEN_USE_MKL_ALL
#include "Field.hpp"
#include "FTLE.hpp"
#include "Helicity.hpp"
#include "hausdorff.hpp"
#include "Interpolation.hpp"
#include <stdlib.h>
#include "Frechet.hpp"

int main(void)
{
    double pi = 4.0 * std::atan(1.0);
    LCS::position<double> domain(600, 600, 1);
    domain.setAll(0., 2.*pi, 0., 2.*pi, -0.012, 0.012);

    domain.setInitialTime(0.); domain.setFinalTime(4.);
    domain.setAbsTol(1e-013);
    domain.setRelTol(1e-013);
    domain.auxiliaryGridSizingFactor = .05;
    // domain.computeFlowProperties();

    LCS::Helicity<double> helicityField(domain);
    std::vector<std::vector<int>> indices;
    std::vector<std::vector<Point<double>>> allStrainlines, filteredStrainlinesFrechet, filteredStrainlinesHausdorff;
    
    // helicityField.computeHelicityField();
    std::cout << "here" << std::endl;
    helicityField.computeHelicityFieldAuxiliaryGrid();

    // for (double threshold = 1e-04; threshold < 1e-03; threshold += 1e-04)
    // {

    //     helicityField.getHelicityBelowThreshold(threshold, indices);
    //     helicityField.integrateStrainLines(indices, threshold); // indices, tolerance
    //     std::string comm = "mkdir -p ../strainlines_" + std::to_string(threshold * 10000) + " && cp -r ../strainlines/ ../strainlines_" + std::to_string(threshold * 10000) + "/";
    //     system(comm.c_str());
    // }

    /* Create reduced field. */
    //LCS::position<double> reducedPosition(300, 5, 1);
    //reducedPosition.setAll(0., 2.*pi, 0., 2.*pi, -0.012, 0.012);
    //reducedPosition.setAbsTol(1e-013);
    //reducedPosition.setRelTol(1e-013);
    //reducedPosition.auxiliaryGridSizingFactor = 0.1;
    //reducedPosition.setInitialTime(0.);
    //reducedPosition.setFinalTime(4.);
    //LCS::Helicity<double> reducedField(reducedPosition);
    //createInterpolatedGrid(helicityField, domain, reducedPosition, reducedField);
    //std::cout << "here " << std::endl;
    //std::cout.flush();
    //reducedField.getHelicityBelowThreshold(2e-04, indices);
    //reducedField.integrateStrainLines(indices, 2e-04, allStrainlines);

    //std::ofstream reducedOutput; reducedOutput.open("../reducedHelicity");
    // std::ofstream output; output.open("../helicity");

    // for (int i = 0; i < domain.getXExtent(); ++i)
    // {
    //     for (int j = 0; j < domain.getYExtent(); ++j)
    //     {
    //         Point<double> thisPoint = domain.getValue(i, j, 0);
    //         output << thisPoint.x << "\t" << thisPoint.y << "\t" << helicityField.helicityField_[i][j][0] << std::endl;
    //     }
    // }

    //for (int i = 0; i < reducedField.getXExtent(); ++i)
    //{
    //    for (int j =  0; j < reducedField.getYExtent(); ++j)
    //    {
    //        Point<double> thisPoint = reducedPosition.getValue(i, j, 0);
    //        reducedOutput << thisPoint.x << "\t" << thisPoint.y << "\t" << reducedField.helicityField_[i][j][0] << std::endl;
    //    }
    //}
    //exit(-1);
    // helicityField.writeHelicityField();
    // reducedField.writeHelicityBelowThreshold(1.e-04);
    helicityField.getHelicityBelowThreshold(1e-04, indices);
    helicityField.integrateStrainLines(indices, 1e-04, allStrainlines); // indices, tolerance
    Frechet::filterOnFrechet(allStrainlines, 0.1, filteredStrainlinesFrechet);
    std::cout << "Frechet found " << filteredStrainlinesFrechet.size() << " filtered strainlines. " << std::endl;
    filterOnHausdorff(filteredStrainlinesFrechet, 0.1, filteredStrainlinesHausdorff);
    Frechet::writeFrechet(filteredStrainlinesHausdorff);
    // reducedField.setXStepIntegration(domain.auxiliaryGridSizingFactor * domain.getXStep());
        // reducedField.setYStepIntegration(domain.auxiliaryGridSizingFactor * domain.getYStep());
    // reducedField.setZStepIntegration(domain.auxiliaryGridSizingFactor * domain.getZStep());
    // helicityField.writeSeedPoints(indices);

    return 0;
}
