#include "observer.hpp"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <ios>
#include <memory>
#include <stdexcept>
#include <string>

aos::simulation::csv_state_observer::csv_state_observer(const std::string& filename, std::size_t num_rods) : m_num_rods(num_rods) {
    // Use Boost.Filesystem to ensure the output directory exists
    boost::filesystem::path file_path(filename);
    if (file_path.has_parent_path()) {
        boost::filesystem::create_directories(file_path.parent_path());
    }

    m_file = std::make_shared<std::ofstream>(filename);
    if (not m_file->is_open()) {
        throw std::runtime_error("Observer could not open output file: " + filename);
    }

    // Set a high precision for floating-point numbers
    const int high_precision = 10;
    *m_file << std::fixed << std::setprecision(high_precision);

    // --- Write the CSV Header ---
    *m_file << "time,q_w,q_x,q_y,q_z,w_x,w_y,w_z";
    for (std::size_t i = 0; i < m_num_rods; ++i) {
        *m_file << ",M_" << (i + 1);
    }
    *m_file << "\n";
}

// NOLINTBEGIN(readability-magic-numbers)
void aos::simulation::csv_state_observer::operator()(const system_state& state, double time) const {
    *m_file << time << ','                         //
            << state.attitude.coeffs().w() << ','  //
            << state.attitude.coeffs().x() << ','  //
            << state.attitude.coeffs().y() << ','  //
            << state.attitude.coeffs().z() << ','  //
            << state.angular_velocity.x() << ','   //
            << state.angular_velocity.y() << ','   //
            << state.angular_velocity.z();

    for (std::size_t i = 0; i < m_num_rods; ++i) {
        *m_file << ',' << state.rod_magnetizations(static_cast<int>(i));
    }
    *m_file << '\n';
}
// NOLINTEND(readability-magic-numbers)
