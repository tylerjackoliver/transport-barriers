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
