#ifndef __FIELD_H__
#define __FIELD_H__

#include <iostream>
#include "Point.hpp"
#include <boost/multi_array.hpp>

template <class Type>
class field
{

    public:

        boost::multi_array<Point<Type>, 3> data;        // 3-D array of points
        
        unsigned xDim;
        unsigned yDim; 
        unsigned zDim;
        
        field(unsigned xDim, unsigned yDim, unsigned zDim);
        ~field();

};

template <typename Type>
field<Type>::field(unsigned _xDim, unsigned _yDim, unsigned _zDim) : xDim(_xDim), yDim(_yDim), zDim(_zDim)
{
    // Initialise data to be the correct size based on xDim/yDim/zDim
    boost::multi_array<Point<Type>, 3>::extent_gen extents;
    this->data.resize(extents[xDim+2][yDim+2][zDim+2]);

    // Set data to be -1-indexed
    this->data.rebase(-1);
}

template <typename Type>
field<Type>::~field(){};

#include "field_new.tpp"

#endif