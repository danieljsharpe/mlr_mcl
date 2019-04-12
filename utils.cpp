#include <random>
#include <typeinfo>

// uniform random number generator
template <class T>
static T rand_unif(T xmin, T xmax, int seed) {
    static std::default_random_engine generator (seed);
    if (typeid(T) == typeid(double)) {
        std::uniform_real_distribution<double> distribution1(xmin,xmax);
        return distribution1(generator); }
    else if (typeid(T) == typeid(int)) {
        std::uniform_int_distribution<int> distribution2(xmin,xmax);
        return distribution2(generator); }
}
