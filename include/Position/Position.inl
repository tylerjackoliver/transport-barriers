/* @brief Set data values of the field using \f$x\f$, \f$y\f$, \f$z\f$ ranges
*
* @param[in] xrange Vector that contains \f$x\f$ coordinates of the grid
* @param[in] yrange Vector that contains \f$y\f$ coordinates of the grid
* @param[in] zrange Vector that contains \f$z\f$ coordinates of the grid
*/
template <typename Type>
void position<Type>::setAll(const std::vector<Type>& xrange, const std::vector<Type>& yrange, const std::vector<Type>& zrange)
{
    /* Make sure the input sizes match */
    if (xrange.size() != (this->nx_ + 2) || yrange.size() != (this->ny_ + 2) || zrange.size() != (this->nz_ + 2) )
    {
        throw std::domain_error("Input ranges in setAll do not match.");
    }

    /* Update the minima and maxima for each direction based on xrange/yrange/zrange */
    this->xMin = xrange[0];
    this->xMax = xrange[xrange.size() - 1];
    this->yMin = yrange[0];
    this->yMax = yrange[yrange.size() - 1];
    this->zMin = zrange[0];
    this->zMax = zrange[zrange.size() - 1];

    /* Add the data from each of these ranges to the main data arrays */
    #pragma omp parallel for /* Maximise throughput */
    for (int i = -1; i < (static_cast<int>(this->nx_) + 1); ++i)
    {
        for (int j = -1; j < static_cast<int>(this->ny_ + 1); ++j)
        {
            for (int k = -1; k < static_cast<int>(this->nz_ + 1); ++k)
            {
                this->data_[i][j][k].x = xrange[i+1];
                this->data_[i][j][k].y = yrange[j+1];
                this->data_[i][j][k].z = zrange[k+1];
            }
        }
    }
}

/* @brief Set data values of the field using \f$x\f$, \f$y\f$, \f$z\f$ coordinates.
*   Manually calls the vector form of setAll during execution.
*
* @param[in] xmin Minimum \f$x\f$ coordinate
* @param[in] xmax Maximum \f$x\f$ coordinate
* @param[in] ymin Minimum \f$y\f$ coordinate
* @param[in] ymax Maximum \f$y\f$ coordinate
* @param[in] zmin Minimum \f$z\f$ coordinate
* @param[in] zmax Maximum \f$z\f$ coordinate
*/
template <typename Type>
void position<Type>::setAll(Type xmin, Type xmax, Type ymin, Type ymax, Type zmin, Type zmax)
{
    /* Manually construct the range vectors to pass into the vector form of setAll */
    std::vector<Type> xrange(this->nx_+2, 0.), yrange(this->ny_+2, 0.), zrange(this->nz_+2, 0.);

    int i = -1;
    int j = -1;
    int k = -1; // due to -1 extent

    if (this->nx_ < 2)
    {
        Type gap = (xmax - xmin) / (this->nx_+1);
        this->xStep = gap;
        i = 0;
        std::generate(xrange.begin(), xrange.end(), \
            [&]{ return xmin + (i++) * gap;});
    } else
    {
        this->xStep = (xmax-xmin) / (this->nx_-1);
        std::generate(xrange.begin(), xrange.end(), \
            [&]{ return xmin + (i++) * (xmax-xmin) / (this->nx_-1);});
    }

    if (this->ny_ < 2)
    {
        Type gap = (ymax - ymin) / (this->ny_+1);
        this->yStep = gap;
        j = 0;
        std::generate(yrange.begin(), yrange.end(), 
            [&]{ return ymin + (j++) * gap;});
    } else
    {
        this->yStep = (ymax-ymin) / (this->ny_-1);
        std::generate(yrange.begin(), yrange.end(), \
            [&]{ return ymin + (j++) * (ymax-ymin) / (this->ny_-1);});
    }

    // Correction for if we're on a single z-plane (planar flow)
    if (this->nz_ < 2)
    {
        Type gap = (zmax - zmin) / (this->nz_+1);
        this->zStep = gap;
        k = 0;
        std::generate(zrange.begin(), zrange.end(), 
            [&]{ return zmin + (k++) * gap;});
    } else
    {
        this->zStep = (zmax-zmin) / (this->nz_-1);
        std::generate(zrange.begin(), zrange.end(), \
            [&]{ return zmin + (k++) * (zmax-zmin) / (this->nz_-1);});
    }
    setAll(xrange, yrange, zrange);
}


/* @brief Advect the initial flow coordinates forward using a Burlisch-Stoer numerical integration scheme
*
* @param absTol Absolute tolerance to use in the integration
* @param relTol Relative tolerance to use in the integration
*/
template <typename Type>
void position<Type>::advectPosition(double absTol, double relTol)
{
    unsigned prog = 0;

    #pragma omp parallel for shared(prog)
    for (int i = -1; i < static_cast<int>(this->nx_+1); i++)
    {   
        for (int j = -1; j < static_cast<int>(this->ny_+1); j++)
        {
            for (int k = -1; k < static_cast<int>(this->nz_+1); k++)
            {
                std::vector<Type> x(3);
                x[0] = this->data_[i][j][k].x;
                x[1] = this->data_[i][j][k].y;
                x[2] = this->data_[i][j][k].z;

                Helpers::integrateForceFunction(x, this->initTime_, this->finalTime_, absTol, relTol);
                
                this->updatedPosition_[i][j][k].x = x[0];
                this->updatedPosition_[i][j][k].y = x[1];
                this->updatedPosition_[i][j][k].z = x[2];
            }
        }
        #pragma omp atomic
        prog++;
        std::cout << "Completed advecting " << double(prog)/(this->nx_+1) * 100. << "% of the flow." << std::endl;
    }
}