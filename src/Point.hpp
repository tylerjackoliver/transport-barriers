#ifndef _POINT_H_
#define _POINT_H_

#include <ostream>
#include <vector>
#include <type_traits>
#include <stdexcept>

/* @brief Custom class that holds coordinates of a given point in the field.
 *
 * This supports an arbitrary number of points using C++17 fold expressions.
 */
template <typename Type, int dimension>
struct Point
{
	/* Underlying storage container. */
	std::vector<Type> state;

	/* Define class constructors */

	/* @brief Default constructor. Leaves everything to be created later but sets the size.
	 */
	Point( void ) : state(dimension);
	{};

	/* @brief Constructor using arbitrary number of parameters.
	 * @param[in] std::initialiser_list<Type> containing the values to assign.
	 */
	Point(std::initialiser_list<Type>& args) : state(args)
    {}

	/* Initialise a point from another point */
    Point(const Point<Type, dimension>& in)
    {
    	this->state = in.state;
    }
    
    /* @brief Return the size (= dimension) of the underlying storage.
     * @returns The size of the underlying storage.
     */
    int size( void )
    {
    	return this->state.size();
    }

    /* Operator overloads */

    Point<Type, dimension>& operator+=(Point<Type, dimension> const& rhs)
    {
    	for (int i = 0; i < dimension; ++i) this->state[i] += rhs.state[i];
    	return *this;
    }

    Point<Type, dimension>& operator-=(Point<Type, dimension> const& rhs)
    {
    	for (int i = 0; i < dimension; ++i) this->state[i] -= rhs.state[i];
        return *this;
    }

    template <typename multType>
    Point<Type, dimension>& operator*=(multType a)
    {
    	for (int i = 0; i < dimension; ++i) this->state[i] *= a;
        return *this;        
    }

    template <typename multType>
    Point<Type, dimension>& operator/=(multType a)
    {
    	for (int i = 0; i < dimension; ++i) this->state[i] /= a;
    	return *this;
    }

    Point<Type, dimension> operator+(Point<Type, dimension> const &a) const
    {
    	Point<Type, dimension> toReturn;
    	for (int i = 0; i < dimension; ++i) toReturn.state[i] = a.state[i] + this->state[i];
    	return toReturn;
    }

    Point<Type, dimension> operator-(Point const &a) const
    {
    	Point<Type, dimension> toReturn;
    	for (int i = 0; i < dimension; ++i) toReturn.state[i] = a.state[i] - this->state[i];
    	return toReturn;
    }

    template <typename multType>
    Point<Type, dimension> operator*(multType a)
    {
    	Point<Type, dimension> toReturn;
    	for (int i = 0; i < dimension; ++i) toReturn.state[i] = a.state[i] * this->state[i];
        return toReturn;
    }

    bool operator==(const Point<Type, dimension>& a) const
    {
        for (int i = 0; i < dimension; ++i)
        {
        	if (a.state[i] != this->state[i]) return false;
        }
        return true;
    }

    void vectorToPoint(const std::vector<Type> &in)
    {
    	if (in.size() != dimension) throw std::runtime_error("The sizes of vector in vectorToPoint do not match.");
    	this->state = in;
    }

};

template <typename Type>
std::ostream& operator<<(std::ostream& os, const Point<Type>& in)
{
    return (os << "(" << in.x << ", " << in.y << ", " << in.z << ")");    
}

#endif // _POINT_H_
