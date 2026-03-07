#pragma once

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
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
    simulation(const std::string&                  output_filename,
               const simulation_properties&        properties,
               const std::shared_ptr<spacecraft>&  satellite,
               const std::shared_ptr<environment>& environment);
    simulation(const simulation_properties& properties,
               std::shared_ptr<spacecraft>  satellite,
               std::shared_ptr<environment> environment,
               std::shared_ptr<dynamics>    dynamics,
               std::shared_ptr<observer>    observer);

    void run();

protected:

    void fix_integration_errors();

private:

    std::shared_ptr<spacecraft>  _satellite;
    std::shared_ptr<environment> _environment;
    std::shared_ptr<dynamics>    _dynamics;
    std::shared_ptr<observer>    _observer;
    system_state                 _current_state;
    real_t                       _t_start;
    real_t                       _t_end;
    real_t                       _t_now;
    real_t                       _dt_initial;
    real_t                       _checkpoint_interval;
    real_t                       _absolute_error;
    real_t                       _relative_error;
    int                          _stepper_function;
};

}  // namespace aos
