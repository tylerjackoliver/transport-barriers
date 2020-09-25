#ifndef _POINT_H_
#define _POINT_H_

#include <ostream>

template <typename Type>
struct Point
{
    Type x;
    Type y;
    Type z;

    // Class constructors

    Point(){};  // Default initialiser;leave x, y, z to be created later
    Point(Type x_, Type y_, Type z_)
    {
        this->x = x_;
        this->y = y_;
        this->z = z_;
    }
    
    Point(const Point& in)
    {
        this->x = in.x;
        this->y = in.y;
        this->z = in.z;
    }

    // Operator overloads

    Point& operator+=(Point const& rhs)
    {
        this->x += rhs.x;
        this->y += rhs.y;
        this->z += rhs.z;

        return *this;
    }

    Point& operator-=(Point const& rhs)
    {
        this->x -= rhs.x;
        this->y -= rhs.y;
        this->z -= rhs.z;

        return *this;
    }

    template <typename multType>
    Point& operator*=(multType a)
    {
        this->x *= a;
        this->y *= a;
        this->z *= a;

        return *this;        
    }

    template <typename multType>
    Point& operator/=(multType a)
    {
        this->x /= a;
        this->y /= a;
        this->z /= a;
    }

    Point operator+(Point const &a) const
    {
        return Point(a.x+this->x, a.y+this->y, a.z+this->z);
    }

    Point operator-(Point const &a) const
    {
        return Point(a.x-this->x, a.y-this->y, a.z-this->z);
    }

    template <typename multType>
    Point operator*(multType a)
    {
        return Point(this->x*a, this->y*a);
    }

    bool operator==(const Point& a) const
    {
        return (this->x == a.x && this->y == a.y && this->z == a.z);
    }

    void vectorToPoint(const std::vector<Type> &in)
    {
        this->x = in[0];
        this->y = in[1];
        this->z = in[2];
    }

};

template <typename Type>
std::ostream& operator<<(std::ostream& os, const Point<Type>& in)
{
    return (os << "(" << in.x << ", " << in.y << ", " << in.z << ")");    
}

#endif // _POINT_H_