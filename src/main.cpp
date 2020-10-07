#include <iostream>
#define EIGEN_USE_MKL_ALL
#include "field.hpp"
#include "FTLE.hpp"
#include "Helicity.hpp"

int main(void)
{
    double pi = 4.0 * std::atan(1.0);
    LCS::position<double> domain(1000, 1000, 1);
    domain.setAll(0., 2.*pi, 0., 2*pi, -.5, .5);

    domain.setInitialTime(0.); domain.setFinalTime(10.);
    domain.setAbsTol(1e-012);
    domain.setRelTol(1e-012);
    domain.useAuxiliaryGrid = false;
    domain.advectPosition(1e-012, 1e-012);
    domain.computeFlowProperties();
    // domain.writeToFile();

    LCS::FTLE<double> ftleField(domain);
    ftleField.useAuxiliaryGrid = false;
    ftleField.auxiliaryGridSizingFactor = 0.25;
    ftleField.computeFTLE();
    // ftleField.writeFTLE();

    // LCS::Helicity<double> helicityField(domain);
    // // helicityField.computeEigenvectorField();
    // helicityField.computeHelicityField();
    // // helicityField.writeEigenVector();
    // helicityField.writeHelicityField();
    // helicityField.writeHelicityBelowThreshold(1.e-02);

    return 0;
}