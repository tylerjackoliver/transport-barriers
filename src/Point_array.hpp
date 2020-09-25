#ifndef _POINT_H_
#define _POINT_H_

#include <iostream>
#include <mkl.h>
#include <vector>

template <class Type>
struct Points
{
    Type *x;    // Data array for 'x' coordinates
    Type *y;    // Data array for 'y' coordinates
    Type *z;    // Data array for 'z' coordinates

    unsigned xDim;
    unsigned yDim;
    unsigned zDim;

    void getCoordinate(unsigned i, unsigned j, unsigned k, Type &coordinate) const;
    void getCoordinate(std::vector<Type> &coordinate) const;

    void allocateMemory();
    void reallocateMemory();

    Points();
    Points(unsigned _xDim, unsigned _yDim, unsigned _zDim);
    ~Points();
    
    // void setCoordinate(Type &coordinate);
    // void setCoordinate(std::vector<Type> &coordinate);
};

#endif