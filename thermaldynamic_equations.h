#ifndef THERMALDYNAMIC_EQUATIONS_H
#define THERMALDYNAMIC_EQUATIONS_H
#include <cmath>

using namespace std;

class thermaldynamic_equations
{
public:
    thermaldynamic_equations();
    static void total_to_static(const double &total_pressure, const double &total_temperature,
                                const double &area, const double &mass_flow,
                                const double &gas_constant, const double &gamma,
                                double &static_pressure, double &static_temperature, double &density,
                                const double &c_t=0.0);
};

#endif // THERMALDYNAMIC_EQUATIONS_H
