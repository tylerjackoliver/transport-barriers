/* Implementation file the Helicity class */


/* @brief Class constructor for the Helicity class.
   @param[in] pos Reference to a populated position class
*/
Helicity::Helicity(position<double>& _pos) : pos_(_pos)
{
    // Correctly initialise the boost::multiarray as public class member
    typedef boost::multi_array<Eigen::Vector3d, 3> tempType;
    typename tempType::extent_gen extents;
    this->eigenvectorField_.resize(extents[pos_.getXExtent()+2][pos_.getYExtent()+2][pos_.getZExtent()+2];
    this->eigenvectorField_.reindex(-1);

}

/* @brief Private member function that determines the sign of an expression
 * @param[in] val Numeric type to determine the sign of
 * @returns integer (+1 or -1) that gives the sign of val
 */
template <typename Type>
int Helicity::sgn(Type val)
{
    return ( (val > 0) - (val < 0) );   
}

/* @brief Computes the helicity field for a given system using finite differencing with or without an auxiliary grid.
*/
void Helicity::computeEigenvectorField()
{
    if (useAuxiliaryGrid) this->_computeHelicityAuxiliaryGrid();
    else this->_computeHelicityFiniteDifferencing();
}

/* @brief Computes the helicity field using an auxiliary grid
*/
template <typename Type>
void Helicity::_computeEigenvectorAuxiliaryGrid()
{
    unsigned int prog = 0;  // Progress indicator
    
#if defined(_OPENMP)
    #pragma omp parallel for shared(prog)
#else
    #pragma omp simd
#endif
    for (unsigned i = 0; i < pos_.getXExtent(); ++i)
    {
        for (unsigned j = 0; j < pos_.getYExtent(); ++j)
        {
            for (unsigned k = 0; k < pos_.getZExtent(); ++k)
            {
                Point<Type> xNext, xPrev, yNext, yPrev, zNext, zPrev;
                Point<Type> x0Next, x0Prev, y0Next, y0PRev, z0Next, z0Prev;
                Eigen::Matrix<Type, 3, 3> deformation, cauchy_green;

                std::vector<double> xPlus(3), xMinus(3), yPlus(3), yMinus(3), zPlus(3), zMinus(3); // std::vector representations of points for integration
                Point<Type> referencePoint(3), previousReferenceX(3), previousReferenceY(3), previousReferenceZ(3);

                // Get the point we're studying and points in each axis around it (boundary conditions OK)
                referencePoint = pos_.getValue(i, j, k);
                previousReferenceX = pos_.getValue(i-1, j, k);
                previousReferenceY = pos_.getValue(i, j-1, k);
                previousReferenceZ = pos_.getValue(i, j, k-1);

                // Get auxiliary grid spacings
                double auxGridSizeX = (referencePoint.x - previousReferenceX.x) * this->auxiliaryGridSizingFactor;
                double auxGridSizeY = (referencePoint.y - previousReferenceY.y) * this->auxiliaryGridSizingFactor;
                double auxGridSizeZ = (referencePoint.z - previousReferenceZ.z) * this->auxiliaryGridSizingFactor;

                // Point -> vector conversion
                reference[0] = referencePoint.x;
                reference[1] = referencePoint.y;
                reference[2] = referencePoint.z;
                
                // Set up grid around reference point
                xPlus = reference; xPlus[0] += auxGridSizeX; x0Next.vectorToPoint(xPlus);
                xMinus = reference; xMinus[0] -= auxGridSizeX; x0Prev.vectorToPoint(xMinus);
                
                yPlus = reference; yPlus[1] += auxGridSizeY; y0Next.vectorToPoint(yPlus);
                yMinus = reference; yMinus[1] -= auxGridSizeY; y0Prev.vectorToPoint(yMinus);

                zPlus = reference; zPlus[2] += auxGridSizeZ; z0Next.vectorToPoint(zPlus);
                zMinus = reference; zMinus[2] -= auxGridSizeZ; z0Prev.vectorToPoint(zMinus);

                // Integrate all of them
                this->integrate(xPlus); xNext.vectorToPoint(xPlus);
                this->integrate(xMinus); xPrev.vectorToPoint(xMinus);
                this->integrate(yPlus); yNext.vectorToPoint(yPlus);
                this->integrate(yMinus); yPrev.vectorToPoint(yMinus);
                this->integrate(zPlus); zNext.vectorToPoint(zPlus);
                this->integrate(zMinus); zPrev.vectorToPoint(zMinus);

                // Deformation tensor
                deformation(0,0) = (xNext.x - xPrev.x) / (x0Next.x - x0Prev.x);
                deformation(0,1) = (yNext.x - yPrev.x) / (y0Next.y - y0Prev.y);
                deformation(0,2) = (zNext.x - zPrev.x) / (z0Next.z - z0Prev.z);

                deformation(1,0) = (xNext.y - xPrev.y) / (x0Next.x - x0Prev.x);
                deformation(1,1) = (yNext.y - yPrev.y) / (y0Next.y - y0Prev.y);
                deformation(1,2) = (zNext.y - zPrev.y) / (z0Next.z - z0Prev.z);

                deformation(2,0) = (xNext.z - xPrev.z) / (x0Next.x - x0Prev.x);
                deformation(2,1) = (yNext.z - yPrev.z) / (y0Next.y - y0Prev.y);
                deformation(2,2) = (zNext.z - zPrev.z) / (z0Next.z - z0Prev.z);

                cauchy_green = deformation.transpose() * deformation;
                this->eigenvectorField[i][j][k] = cauchy_green.template selfadjointView<Eigen::Lower>().eigenvectors.col(2).real(); // Eigenvalues sorted in ascending order, so want col(2) for biggest eigenvalues
                
                } // k

            } // j

#if defined(_OPENMP)
            #pragma omp atomic
            prog++;
#else
            prog++;
#endif
       } // i

}; // function


/* @brief Computes the helicity field using finite differencing
 */
template <typename Type>
void Helicity::_computeEigenvectorAuxiliaryGrid()
{
    unsigned int prog = 0;  // Progress indicator

#if defined(_OPENMP)
#pragma omp parallel for shared(prog)
#else
#pragma omp simd
#endif
    for (unsigned i = 0; i < pos_.getXExtent(); ++i)
    {
        for (unsigned j = 0; j < pos_.getYExtent(); ++j)
        {
            for (unsigned k = 0; k < pos_.getZExtent(); ++k)
            {
                Point<Type> xNext, xPrev, yNext, yPrev, zNext, zPrev;
                Point<Type> x0Next, x0Prev, y0Next, y0PRev, z0Next, z0Prev;
                Eigen::Matrix<Type, 3, 3> deformation, cauchy_green;

                std::vector<double> xPlus(3), xMinus(3), yPlus(3), yMinus(3), zPlus(3), zMinus(3); // std::vector representations of points for integration
                Point<Type> referencePoint(3), previousReferenceX(3), previousReferenceY(3), previousReferenceZ(3);


                // Set up grid around reference point
    
                // Integrate all of them
                this->integrate(xPlus); xNext.vectorToPoint(xPlus);
                this->integrate(xMinus); xPrev.vectorToPoint(xMinus);
                this->integrate(yPlus); yNext.vectorToPoint(yPlus);
                this->integrate(yMinus); yPrev.vectorToPoint(yMinus);
                this->integrate(zPlus); zNext.vectorToPoint(zPlus);
                this->integrate(zMinus); zPrev.vectorToPoint(zMinus);

                // Deformation tensor
                deformation(0,0) = (xNext.x - xPrev.x) / (x0Next.x - x0Prev.x);
                deformation(0,1) = (yNext.x - yPrev.x) / (y0Next.y - y0Prev.y);
                deformation(0,2) = (zNext.x - zPrev.x) / (z0Next.z - z0Prev.z);

                deformation(1,0) = (xNext.y - xPrev.y) / (x0Next.x - x0Prev.x);
                deformation(1,1) = (yNext.y - yPrev.y) / (y0Next.y - y0Prev.y);
                deformation(1,2) = (zNext.y - zPrev.y) / (z0Next.z - z0Prev.z);

                deformation(2,0) = (xNext.z - xPrev.z) / (x0Next.x - x0Prev.x);
                deformation(2,1) = (yNext.z - yPrev.z) / (y0Next.y - y0Prev.y);
                deformation(2,2) = (zNext.z - zPrev.z) / (z0Next.z - z0Prev.z);

                cauchy_green = deformation.transpose() * deformation;
                this->eigenvectorField[i][j][k] = cauchy_green.template selfadjointView<Eigen::Lower>().eigenvectors().col(2).real(); // Eigenvalues sorted in ascending order, so want col(2) for biggest eigenvalues

            } // k

        } // j

        prog++;
#if defined(_OPENMP)
#pragma omp atomic
        prog++;
#else
        prog++;
#endif
    } // i

}; // function


/* @brief Computes the helicity field using an auxiliary grid
 */
template <typename Type>
void Helicity::_computeEigenvectorFiniteDifferencing()
{
    unsigned prog = 0;
    #pragma omp parallel for shared(prog)
    for (unsigned i = 0; i < pos_.getXExtent(); ++i)
    {
        Point<Type> xNext, xPrev, yNext, yPrev, zNext, zPrev;
        Point<Type> x0Next, x0Prev, y0Next, y0Prev, z0Next, z0Prev;
        Eigen::Matrix<Type, 3, 3> deformation, cauchy_green;
        for (unsigned j = 0; j < pos_.getYExtent(); ++j)
        {
            for (unsigned k = 0; k < pos_.getZExtent(); ++k){
                
                xPrev = pos_.getUpdatedValue(i-1, j, k);
                xNext = pos_.getUpdatedValue(i+1, j, k);
                
                yPrev = pos_.getUpdatedValue(i,j-1,k);
                yNext = pos_.getUpdatedValue(i, j+1, k);
                
                zPrev = pos_.getUpdatedValue(i, j, k-1);
                zNext = pos_.getUpdatedValue(i, j, k+1);
                
                x0Prev = pos_.getValue(i-1, j, k);
                x0Next = pos_.getValue(i+1, j, k);
                
                y0Prev = pos_.getValue(i, j-1, k);
                y0Next = pos_.getValue(i, j+1, k);
                
                z0Prev = pos_.getValue(i, j, k-1);
                z0Next = pos_.getValue(i, j, k+1);

                // deformation tensor 
                deformation(0,0) = (xNext.x-xPrev.x) / (x0Next.x-x0Prev.x);
                deformation(0,1) = (yNext.x-yPrev.x) / (y0Next.y-y0Prev.y);
                deformation(0,2) = (zNext.x-zPrev.x) / (z0Next.z-z0Prev.z);
                
                deformation(1,0) = (xNext.y-xPrev.y) / (x0Next.x-x0Prev.x);
                deformation(1,1) = (yNext.y-yPrev.y) / (y0Next.y-y0Prev.y);
                deformation(1,2) = (zNext.y-zPrev.y) / (z0Next.z-z0Prev.z);

                deformation(2,0) = (xNext.z-xPrev.z) / (x0Next.x-x0Prev.x);
                deformation(2,1) = (yNext.z-yPrev.z) / (y0Next.y-y0Prev.y);
                deformation(2,2) = (zNext.z-zPrev.z) / (z0Next.z-z0Prev.z);

                cauchy_green = deformation.transpose() * deformation;
                auto eivals = 
                eigenvectorField_[i][j][k] = cauchy_green.template selfadjointView<Eigen::Lower>().eigenvectors().col(2).real();

            }
            
        }
    
        prog++;
        std::cout << "Completed processing " << double(prog)/this->nx_ * 100. << "% of the FTLE." << std::endl;
    }
}



