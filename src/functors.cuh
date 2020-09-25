
struct abcSystem{

    struct abcFunctor
    {

        template<class T>
        __host__ __device__ 
        void operator()(T t) const{

            // Unpack the parameters
            value_type x = thrust::get<0>(t);
            value_type y = thrust::get<1>(t);
            value_type z = thrust::get<2>(t);

            thrust::get<3>(t) = sqrt(3.0) * sin(z) + 1.0 * cos(y);
            thrust::get<4>(t) = sqrt(2.0) * sin(x) + sqrt(3.0) * cos(z);
            thrust::get<5>(t) = 1.0 * sin(y) + sqrt(2) * cos(x);

        }

    };

    abcSystem(size_t N) : dim_N(N){}

    template <class state, class deriv>
    void operator()(const state& x, deriv& xdot, value_type t) const
    {

        thrust::for_each(

            thrust::make_zip_iterator( thrust::make_tuple(

                boost::begin(x),
                boost::begin(x) + dim_N,
                boost::begin(x) + 2 * dim_N,
                boost::begin(xdot),
                boost::begin(xdot) + dim_N,
                boost::begin(xdot) + 2 * dim_N)
                
            ),

            thrust::make_zip_iterator(thrust::make_tuple(

                boost::begin(x) + dim_N,
                boost::begin(x) + 2 * dim_N,
                boost::begin(x) + 3 * dim_N,
                boost::begin(xdot) + dim_N,
                boost::begin(xdot) + 2 * dim_N,
                boost::begin(xdot) + 3 * dim_N )

            ), abcFunctor() );

    };

    size_t dim_N;

};