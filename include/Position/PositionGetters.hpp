
/* @brief Getter function for the eigenvalue field
    @param[in] i x-coordinate of the desired point
    @param[in] j y-coordinate of the desired point
    @param[in] k z-coordinate of the desired point
    @returns The eigenvalue at that point
*/
template <typename Type>
double position<Type>::getEigenvalue(int i, int j, int k)
{
    return this->eigenvalueField_[i][j][k];
}

/* @brief Getter function for the eigenvector field
    @param[in] i x-coordinate of the desired point
    @param[in] j y-coordinate of the desired point
    @param[in] k z-coordinate of the desired point
    @returns The eigenvector at that point
*/
template <typename Type>
Eigen::Vector3d position<Type>::getEigenvector(int i, int j, int k)
{
    return this->eigenvectorField_[i][j][k];
}

/* @brief Getter function for the x-step of the position field
   @returns The x-step of the position field.
*/
template <typename Type>
double position<Type>::getXStep() const
{
    return this->xStep;
}

/* @brief Getter function for the y-step of the position field
   @returns The y-step of the position field.
*/
template <typename Type>
double position<Type>::getYStep() const
{
    return this->yStep;
}

/* @brief Getter function for the z-step of the position field
   @returns The z-step of the position field.
*/
template <typename Type>
double position<Type>::getZStep() const
{
    return this->zStep;
}

/* @brief Get the number of points in the 'x'-direction
*
* @returns The number of points in the 'x'-direction
*/
template <typename Type>
auto position<Type>::getXExtent() const
{
    return this->nx_;
}

/* @brief Get the number of points in the 'y'-direction
*
* @returns The number of points in the 'y'-direction
*/
template <typename Type>
auto position<Type>::getYExtent() const
{
    return this->ny_;
}

/* @brief Get the number of points in the 'z'-direction
*
* @returns The number of points in the 'z'-direction
*/
template <typename Type>
auto position<Type>::getZExtent() const
{
    return this->nz_;
}