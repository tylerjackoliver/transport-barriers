// #include "field.hpp"
#include "cudaTest.cuh"
// #include "ftle.hpp"

void scheibe(void);

int main(void){

    LCS::position domain(251, 251, 251); // Initialise domain with 250 in each direction

    // // Check the middle of the domain

    domain.setAll(0., 2., 0., 1., 0.0, 0.5); // 0-> 2 in every dimension

    // xyz middle = domain.getValue(125, 125, 125);

    // std::cout << "Middle of the domain is " << middle.x << std::endl;

    domain.setInitialTime(0.); domain.setFinalTime(10.);

    domain.advectPosition(1e-08, 1e-08);

    domain.writeToFile();

    // domain.loadFromFile("DGintegration.data2");

    LCS::FTLE ftleField(domain);

    LCS::helicity helo(domain);

    helo.getDominantVectorField();
    helo.computeHelicity();
    helo.writeHelicityField();
    helo.writeHelicityBelowThreshold(1.e-4);
    std::vector<xyz> helicityIndices = helo.getHelicityBelowThreshold(1.e-04);
    helo.propagateStrainlines(helicityIndices);

    // Compute ftle

    // ftleField.computeFTLE();
    // ftleField.writeFTLE();
}