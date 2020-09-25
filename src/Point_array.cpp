#include "Point_array.hpp"

/*

Data array allocation using mkl_calloc. Depending on the architecture being used, we align our
data differently to ensure maximum performance when/if using SSE/AVX+ instructions.

The final argument to mkl_calloc() is the alignment size, *but* ICPC may also choose to raise
the alignment if required.

*/
template <typename Type>
Points<Type>::Points(unsigned _xDim, unsigned _yDim, unsigned _zDim) : xDim(_xDim), yDim(_yDim), zDim(_zDim)
{

    // Allocate data arrays, plus two for boundary conditions
    x = (Type *) mkl_calloc( (xDim + 2), sizeof(Type), 32);
    y = (Type *) mkl_calloc( (yDim + 2), sizeof(Type), 32);
    z = (Type *) mkl_calloc( (zDim + 2), sizeof(Type), 32);

};

template <typename Type>
Points<Type>::Points()
{
    x = NULL;
    y = NULL;
    z = NULL;

    xDim = -1;
    yDim = -1;
    zDim = -1;
}

template <typename Type>
Points<Type>::~Points()
{
    mkl_free(this->x);
    mkl_free(this->y);
    mkl_free(this->z);
}


template <typename Type>
void Points<Type>::allocateMemory()
{

    // Check we're about to allocate a finite size
    if (xDim < 0 || yDim < 0 || zDim < 0)
    {
        std::cout << "Error: one of the array lengths has not been allocated." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Allocate memory
    x = (Type *)mkl_calloc( (xDim + 2), sizeof(Type), 32);
    y = (Type *)mkl_calloc( (yDim + 2), sizeof(Type), 32);
    z = (Type *)mkl_calloc( (zDim + 2), sizeof(Type), 32);

}

template <typename Type>
void Points<Type>::reallocateMemory()
{

    // Remove existing memory
    mkl_free(x);
    mkl_free(y);
    mkl_free(z);

    allocateMemory(); // Re-allocate memory

}