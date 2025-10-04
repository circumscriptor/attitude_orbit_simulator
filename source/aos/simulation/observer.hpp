#pragma once

#include "aos/core/state.hpp"

#include <cstddef>
#include <fstream>
#include <memory>
#include <string>

namespace aos::simulation {

// Use established type aliases for consistency
using core::system_state;

/**
 * @class csv_state_observer
 * @brief An observer functor that writes the simulation state to a CSV file.
 *
 * This class is designed to be called by a boost::odeint integrator at each
 * successful time step. It handles file creation, directory management,
 * header writing, and formatted data logging.
 */
class csv_state_observer {
public:

    /**
     * @brief Constructs the CSV observer.
     * @param filename The path to the output CSV file.
     * @param num_rods The number of hysteresis rods in the simulation, used to
     *                 dynamically generate the header.
     */
    explicit csv_state_observer(const std::string& filename, std::size_t num_rods);

    /**
     * @brief The main operator called by the ODE solver at each step.
     * @param state The current state vector (x) of the system.
     * @param time The current simulation time in seconds.
     */
    void operator()(const system_state& state, double time) const;

private:

    std::shared_ptr<std::ofstream> m_file;
    size_t                         m_num_rods;
};

}  // namespace aos::simulation
