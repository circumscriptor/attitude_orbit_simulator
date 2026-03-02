#include "verification_observer.hpp"

#include "aos/components/spacecraft.hpp"
#include "aos/environment/environment.hpp"
#include "aos/simulation/observer.hpp"
#include "aos/verify/details/verification_observer_impl.hpp"

#include <memory>
#include <string>

namespace aos {

auto verification_observer::create(const std::string&                 filename,
                                   std::shared_ptr<const spacecraft>  satellite,
                                   std::shared_ptr<const environment> environment,
                                   const observer_properties&         properties) -> std::shared_ptr<observer> {
    return std::make_shared<verification_observer_impl>(filename, satellite, environment, properties);
}

}  // namespace aos
