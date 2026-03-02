#include "aos/components/spacecraft.hpp"
#include "aos/environment/environment.hpp"
#include "aos/simulation/observer.hpp"

#include <memory>
#include <string>

namespace aos {

class verification_observer {
public:

    static auto create(const std::string&                 filename,
                       std::shared_ptr<const spacecraft>  satellite,
                       std::shared_ptr<const environment> environment,
                       const observer_properties&         properties) -> std::shared_ptr<observer>;
};

}  // namespace aos
