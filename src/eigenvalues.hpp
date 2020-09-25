template <typename T>
void position::getDominantEigenvalue(Points<T> &data)
{
    MKL_INT matrixLayout = LAPACK_ROW_MAJOR;
    char jobz = 'N';                                                                    // Want only eigenvalues
    char range = 'I';                                                                   // Want only some eigenvalues
    char uplo = 'U';                                                                    // Store the upper triangular portion of A
    MKL_INT n = 3;                                                                      // Order of the matrix
    MKL_INT lda = 3;                                                                    // Size of upper triangular
    double vl = 0.;                                                                     // Not used.
    double vu = 0.;                                                                     // Not used.
    MKL_INT il = 2;                                                                     // Want 2nd-indexed eigenvalue (max here)
    MKL_INT iu = 3;                                                                     // Want 2nd-indexed eigenvalue (max here)
    double absTol = LAPACKE_dlamch('S');                                                // Want accuracy to safe minimum such that 1/S does not overflow (> eps)
    
    MKL_INT ldz = 3;
    MKL_INT info;
    MKL_INT m;                                                                          // Number of eigenvalues
    MKL_INT isuppz[2*n];                                                                // Not used.

    // Allocate memory required
    double *a = (double *) mkl_calloc(lda * n, sizeof(double), 64);                     // Matrix
    double *w = (double *) mkl_calloc(3, sizeof(double), 64);                           // Eigenvalues
    double *z = (double *) mkl_calloc(ldz * (iu - il) * matrixSize, sizeof(double), 64); // Select eigenvalues

    info = LAPACKE_dsyevr( LAPACK_COL_MAJOR, 'V', 'I', 'U', n, a, lda, vl, vu, il, iu, absTol, &m, w, z, ldz, isuppz);

    if (info > 0)
    {
        printf("Convergence failed.\n");
        exit(EXIT_FAILURE);
    }

    eigenvalue = z[0];

}