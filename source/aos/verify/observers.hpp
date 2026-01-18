#pragma once

#include "aos/core/state.hpp"
#include "aos/verify/hysteresis_loop_dynamics.hpp"

#include <fstream>
#include <memory>
#include <string>

namespace aos {

// Observer to write H, M, and B data to a file
struct bh_observer {
public:

    explicit bh_observer(const std::string& filename);
    void operator()(const hysteresis_state_type& m, double t) const;

private:

    std::shared_ptr<std::ofstream> _file;
};

class orbit_observer {
public:

    explicit orbit_observer(const std::string& filename);
    void operator()(const system_state& state, double t) const;

private:

    std::shared_ptr<std::ofstream> _file;
};

class attitude_observer {
public:

    explicit attitude_observer(const std::string& filename);
    void operator()(const system_state& state, double t) const;

protected:

    static auto calculate_nadir_error(const system_state& state) -> double;

private:

    std::shared_ptr<std::ofstream> _file;
};

}  // namespace aos
