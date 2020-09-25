struct xyz
{

    public:

        double x;
        double y;
        double z;

    template <typename> type
    double operator(type i)
    {

        if (i == 0) return this->x;
        if (i == 1) return this->y;
        if (i == 2) return this->z;

    }
 
};
