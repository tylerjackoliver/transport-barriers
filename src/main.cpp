#include <iostream>
// #define EIGEN_USE_MKL_ALL
#include "Field.hpp"
#include "FTLE.hpp"
#include "Helicity.hpp"

int main(void)
{
    double pi = 4.0 * std::atan(1.0);
    LCS::position<double> domain(500, 500, 1);
    domain.setAll(0., 2.*pi, 0., 2.*pi, -0.012, 0.012);

    domain.setInitialTime(0.); domain.setFinalTime(3.);
    domain.setAbsTol(1e-013);
    domain.setRelTol(1e-013);
    domain.auxiliaryGridSizingFactor = 0.1;
    domain.computeFlowProperties();

    LCS::Helicity<double> helicityField(domain);
    std::vector<std::vector<int>> indices;
    
    helicityField.computeHelicityField();
    helicityField.writeHelicityField();
    helicityField.writeHelicityBelowThreshold(1.e-04);
    helicityField.getHelicityBelowThreshold(1e-04, indices);
    helicityField.writeSeedPoints(indices);
    helicityField.filterSeedPoints(indices);
    helicityField.integrateStrainLines(indices, 1e-04); // indices, tolerance

    return 0;
}