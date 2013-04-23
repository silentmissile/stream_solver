#include "thermaldynamic_equations.h"

thermaldynamic_equations::thermaldynamic_equations()
{
}

void thermaldynamic_equations::total_to_static(const double &total_pressure, const double &total_temperature,
                                               const double &area, const double &mass_flow,
                                               const double &gas_constant, const double &gamma,
                                               double &static_pressure, double &static_temperature,
                                               double &density, const double &c_t)
{
    double residence, Cp=gamma*gas_constant/(gamma-1), flow_speed, tmp;
    static_pressure=total_pressure;
    static_temperature=total_temperature;
    do
    {
        tmp=static_pressure;
        density=static_pressure/(gas_constant*static_temperature);
        flow_speed=mass_flow/(density*area);
        static_pressure=total_pressure-2*(std::pow(flow_speed,2)+std::pow(c_t,2));
        static_temperature=total_temperature-0.5*(std::pow(flow_speed,2)+std::pow(c_t,2))/Cp;
        residence=(static_temperature-tmp)/tmp;
    }while(std::fabs(residence)>0.001);
    density=static_pressure/(gas_constant*static_temperature);
}
