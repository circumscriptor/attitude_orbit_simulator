#pragma once

#include "aos/verify/hysteresis_loop_dynamics.hpp"

#include <fstream>
#include <memory>
#include <string>

namespace aos::verify {

// Observer to write H, M, and B data to a file
struct bh_observer {
public:

    explicit bh_observer(const std::string& filename);

    void operator()(const hysteresis_state_type& m, double t) const;

private:

    std::shared_ptr<std::ofstream> _file;
};

}  // namespace aos::verify
