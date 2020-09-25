#include "Point_array.cpp"
#include <iostream>
#include <mkl.h>
#include "domain.hpp"

int main (void)
{
    domain test_domain(251, 251, 251);
    test_domain.initialisePosition(-1, 1, -1, 1, -1, 1);
    std::cout << test_domain.initialDomain.x[0] << std::endl;
}