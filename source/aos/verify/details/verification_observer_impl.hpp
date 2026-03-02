#pragma once

#include "aos/components/spacecraft.hpp"
#include "aos/core/state.hpp"
#include "aos/environment/environment.hpp"
#include "aos/simulation/details/observer_impl.hpp"
#include "aos/simulation/observer.hpp"

#include <memory>
#include <ostream>
#include <string>

namespace aos {

class verification_observer_impl : public observer_impl {
public:

    verification_observer_impl(const std::string&                 filename,
                               std::shared_ptr<const spacecraft>  sat,
                               std::shared_ptr<const environment> env,
                               const observer_properties&         props);

    auto write_header() -> std::ostream& override;
    auto write(const system_state& state, double time) -> std::ostream& override;

private:

    std::shared_ptr<const spacecraft>  _sat;
    std::shared_ptr<const environment> _env;
};

}  // namespace aos
