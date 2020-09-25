#ifndef __EIGENVECTORS_H__
#define __EIGENVECTORS_H__

template <class Type>
class eigenvectorField
{
    public:

        eigenvectorField(Points _domain);
        eigenvectorField(unsigned xDim, unsigned yDim, unsigned zDim);
        ~eigenvectorField();

        void computeEigenvectorField();
        void writeEigenvectorField();
        
        Points<Type> data;
        
    private:

}

#endif
