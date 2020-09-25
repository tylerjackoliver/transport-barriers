#pragma once
// struct xyz
// {

//     double x;
//     double y;
//     double z;

// };

// Type definitions

typedef std::vector<double> state_type;
// Define the integrator type
typedef boost::multi_array<xyz, 3> field_;
typedef boost::multi_array<double, 3> doubleField_;
typedef Eigen::Vector3d eigenvectorType;
typedef boost::multi_array<eigenvectorType, 3> eigenvectorField_;
typedef boost::multi_array_types::extent_range range;
field_::extent_gen extents;
eigenvectorField_::extent_gen evExtents;

// // CUDA definitions
// typedef double value_type;
// typedef thrust::device_vector<value_type> state_typeGPU;
// void advectPositionGPUDriver(field_& initPos, field_& finalPos); // CUDA integrator
