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

    explicit hysteresis_observer(const std::string& filename, real h_max, real frequency);
    void operator()(const hysteresis_state_type& m, real t) const;

private:

    std::shared_ptr<std::ofstream> _file;
    real                           _h_max;
    real                           _frequency;
};

}  // namespace aos
