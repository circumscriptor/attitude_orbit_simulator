#pragma once

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/environment/environment.hpp"
#include "aos/simulation/config.hpp"
#include "aos/simulation/dynamics.hpp"
#include "aos/simulation/observer.hpp"

#include <memory>
#include <string>

namespace aos {

class simulation {
public:

    simulation(const std::string& output_filename, const simulation_properties& properties);
    simulation(const simulation_properties& properties,
               std::shared_ptr<spacecraft>  sat,
               std::shared_ptr<environment> env,
               std::shared_ptr<dynamics>    dyn,
               std::shared_ptr<observer>    obs);

    void run();

protected:

    void fix_integration_errors();

private:

    std::shared_ptr<spacecraft>  _satellite;
    std::shared_ptr<environment> _environment;
    std::shared_ptr<dynamics>    _dynamics;
    std::shared_ptr<observer>    _observer;
    system_state                 _current_state;
    double                       _t_start;
    double                       _t_end;
    double                       _t_now;
    double                       _dt_initial;
    double                       _checkpoint_interval;
    double                       _absolute_error;
    double                       _relative_error;
    int                          _stepper_function;
};

}  // namespace aos
