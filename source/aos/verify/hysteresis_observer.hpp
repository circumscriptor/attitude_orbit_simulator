#pragma once

#include "aos/core/types.hpp"
#include "aos/verify/hysteresis_loop_dynamics.hpp"

#include <fstream>
#include <memory>
#include <string>

namespace aos {

// Observer to write H, M, and B data to a file
struct hysteresis_observer {
public:

    explicit hysteresis_observer(const std::string& filename, real_t h_max, real_t frequency);
    void operator()(const hysteresis_state_type& m, real_t t) const;

private:

    std::shared_ptr<std::ofstream> _file;
    real_t                         _h_max;
    real_t                         _frequency;
};

}  // namespace aos
