#include <Eigen/Dense>
#include <GeographicLib/MagneticModel.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/constants.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace aos {

using mat3x3               = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
using rod_magnetizations_t = Eigen::VectorX<double>;
using quat                 = Eigen::Quaternion<double>;
using vec3                 = Eigen::Matrix<double, 3, 1>;
using inertia_t            = mat3x3;
using velocity_t           = vec3;

constexpr double vacuum_permeability = 4.0 * 3.14159265358979323846e-7;

class permanent_magnet {
public:

    explicit permanent_magnet(vec3 moment_body) : m_magnetic_moment_body(std::move(moment_body)) {}

    [[nodiscard]]
    vec3 get_magnetic_moment() const {
        return m_magnetic_moment_body;
    }

private:

    vec3 m_magnetic_moment_body;  // A·m^2
};

class hysteresis_rod {
public:

    hysteresis_rod(double vol, const vec3& orientation) : m_volume(vol), m_orientation_body(orientation.normalized()) {}

    [[nodiscard]]
    vec3 get_magnetic_moment(double magnetization_state) const {
        // M_vec = (M_scalar / μ₀) * V * b_unit
        return (magnetization_state / vacuum_permeability) * m_volume * m_orientation_body;
    }

    // THIS IS THE PLACEHOLDER for a real hysteresis model (e.g., Jiles-Atherton)
    // It calculates the rate of change of the rod's scalar magnetization
    [[nodiscard]]
    double calculate_magnetization_derivative(double current_magnetization, const vec3& b_field_body) const {
        // H = B ⋅ b (component of B along rod axis)
        double h_along_rod = b_field_body.dot(m_orientation_body);

        double damping_coeff = 0.1;  // Fictitious coefficient
        return damping_coeff * (h_along_rod - current_magnetization);
    }

    [[nodiscard]]
    vec3 get_orientation() const {
        return m_orientation_body;
    }

private:

    double m_volume;            // m^3
    vec3   m_orientation_body;  // Unit vector in body frame
};

class spacecraft {
public:

    explicit spacecraft(const mat3x3& inertia) : m_inertia_tensor(inertia), m_inertia_tensor_inverse(inertia.inverse()) {}

    void add_magnet(const permanent_magnet& magnet) { m_magnets.push_back(magnet); }
    void add_rod(const hysteresis_rod& rod) { m_rods.push_back(rod); }

    // Public accessors
    [[nodiscard]] const mat3x3&                        get_inertia_tensor() const { return m_inertia_tensor; }
    [[nodiscard]] const mat3x3&                        get_inertia_tensor_inverse() const { return m_inertia_tensor_inverse; }
    [[nodiscard]] const std::vector<permanent_magnet>& get_magnets() const { return m_magnets; }
    [[nodiscard]] const std::vector<hysteresis_rod>&   get_rods() const { return m_rods; }
    [[nodiscard]] size_t                               get_num_rods() const { return m_rods.size(); }

private:

    mat3x3                        m_inertia_tensor;
    mat3x3                        m_inertia_tensor_inverse;
    std::vector<permanent_magnet> m_magnets;
    std::vector<hysteresis_rod>   m_rods;
};

struct simulation_parameters {
    inertia_t inertia_tensor;
    vec3      initial_angular_velocity;
    vec3      permanent_magnet_moment;

    double t_start{};
    double t_end{};
    double dt_initial{};

    double orbit_altitude_km{};
    double orbit_inclination_deg{};

    std::string output_filename;
};

struct system_state {
    quat                 attitude;
    velocity_t           angular_velocity;
    rod_magnetizations_t rod_magnetizations;

    explicit system_state(int num_rods = 0) : attitude(quat::Identity()), angular_velocity(vec3::Zero()) {
        if (num_rods > 0) {
            rod_magnetizations.resize(num_rods);
            rod_magnetizations.setZero();
        }
    }
};

class spacecraft_dynamics {
public:

    explicit spacecraft_dynamics(const spacecraft& sc) : m_spacecraft(sc) {}

    // The main function for the ODE solver
    // It computes dxdt = f(x, t)
    void operator()(const system_state& x, system_state& dxdt, double t) const {
        // --- 1. Get External Input: Earth's Magnetic Field ---
        // In a real simulation, this would come from an orbit propagator and a field model.
        // For this example, let's use a simple, rotating field.
        vec3 b_inertial;
        b_inertial << 2e-5 * cos(0.001 * t), 2e-5 * sin(0.001 * t), -4e-5;

        // --- 2. Attitude Kinematics (dq/dt) ---
        const auto& q     = x.attitude;
        const auto& omega = x.angular_velocity;

        quat omega_q;
        omega_q.w()   = 0;
        omega_q.vec() = omega;

        quat q_dot = q * omega_q;
        q_dot.coeffs() *= 0.5;
        dxdt.attitude = q_dot;

        // --- 3. Rotational Dynamics (dω/dt) ---
        // Transform B-field from inertial to body frame
        mat3x3 r_inertial_to_body = q.toRotationMatrix().transpose();
        vec3   b_body             = r_inertial_to_body * b_inertial;

        // Calculate total torque
        vec3 total_torque_body = vec3::Zero();

        // a) Permanent magnet torque
        for (const auto& magnet : m_spacecraft.get_magnets()) {
            total_torque_body += magnet.get_magnetic_moment().cross(b_body);
        }

        // b) Hysteresis rod torque
        for (size_t i = 0; i < m_spacecraft.get_num_rods(); ++i) {
            const auto& rod        = m_spacecraft.get_rods()[i];
            double      m_i        = x.rod_magnetizations(static_cast<int>(i));
            vec3        rod_moment = rod.get_magnetic_moment(m_i);
            total_torque_body += rod_moment.cross(b_body);
        }

        // c) Gyroscopic torque
        total_torque_body -= omega.cross(m_spacecraft.get_inertia_tensor() * omega);

        // Final angular acceleration
        dxdt.angular_velocity = m_spacecraft.get_inertia_tensor_inverse() * total_torque_body;

        // --- 4. Hysteresis Rod Dynamics (dM/dt) ---
        dxdt.rod_magnetizations.resize(m_spacecraft.get_num_rods());
        for (size_t i = 0; i < m_spacecraft.get_num_rods(); ++i) {
            const auto& rod                              = m_spacecraft.get_rods()[i];
            double      m_i                              = x.rod_magnetizations(static_cast<int>(i));
            dxdt.rod_magnetizations(static_cast<int>(i)) = rod.calculate_magnetization_derivative(m_i, b_body);
        }
    }

private:

    const spacecraft& m_spacecraft;
};

struct streaming_observer {
    std::ofstream* output;

    explicit streaming_observer(std::ofstream& output) : output(&output) {}

    void operator()(const system_state& state, double time) const {
        (*output) << time << ',';
        (*output) << state.attitude << ',';
        (*output) << state.angular_velocity << ',';
        (*output) << state.rod_magnetizations << ',';
        (*output) << "\n";
    }
};

static void solve(const spacecraft& satellite, const simulation_parameters& params) {
    using boost::numeric::odeint::integrate_adaptive;
    using boost::numeric::odeint::make_controlled;

    try {
        GeographicLib::MagneticModel mag("wmm2020");
        spacecraft_dynamics          dynamics(satellite /*, mag*/);

        // 3. Set Initial Conditions
        system_state initial_state(static_cast<int>(satellite.get_num_rods()));
        initial_state.angular_velocity = params.initial_angular_velocity;  // Initial tumble in rad/s

        // 4. Setup Simulation and Logging
        std::ofstream      data_file("simulation_output.csv");
        streaming_observer observer(data_file);

        // 5. Run the Integrator
        std::cout << "Starting simulation..." << '\n';

        // Use an adaptive step-size Runge-Kutta-Dormand-Prince 5 method
        using stepper_type = boost::numeric::odeint::runge_kutta_dopri5<system_state>;
        double t_start     = params.t_start;
        double t_end       = params.t_end;  // Simulate for 3000 seconds
        double dt_initial  = params.dt_initial;

        integrate_adaptive(make_controlled<stepper_type>(1e-6, 1e-6), dynamics, initial_state, t_start, t_end, dt_initial, observer);

        std::cout << "Simulation finished. Data saved to simulation_output.csv" << '\n';

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        std::cerr << "Please ensure you have installed the WMM2020 data with: sudo geographiclib-get wmm2020\n";
        return;
    }
}

void generate_default_config(const std::string& config_file_path) {
    std::ofstream ofs(config_file_path);
    if (!ofs) {
        std::cerr << "Error: Could not write to default config file: " << config_file_path << '\n';
        return;
    }

    ofs << "# Configuration for the spacecraft simulation\n"
        << "# All values must be space-separated where multiple tokens are expected.\n\n"

        << "# Inertia tensor [kg m^2] (9 space-separated values, row-major)\n"
        << "inertia_tensor = 0.1 0.0 0.0  0.0 0.2 0.0  0.0 0.0 0.15\n\n"

        << "# Initial angular velocity [rad/s] (3 space-separated values for x, y, z)\n"
        << "initial_angular_velocity = 1.0 -0.5 0.8\n\n"

        << "# Permanent magnet moment [A m^2] (3 space-separated values for x, y, z)\n"
        << "permanent_magnet_moment = 0.0 0.0 1.0\n\n"

        << "# Simulation start time [s]\n"
        << "simulation_t_start = 0.0\n\n"

        << "# Simulation end time [s]\n"
        << "simulation_t_end = 3000.0\n\n"

        << "# Initial time step for adaptive integrator [s]\n"
        << "simulation_dt_initial = 0.1\n\n"

        << "# Circular orbit altitude [km]\n"
        << "orbit_altitude_km = 500.0\n\n"

        << "# Orbit inclination [deg]\n"
        << "orbit_inclination_deg = 51.6\n\n"

        << "# Output CSV filename\n"
        << "output_filename = output/simulation_run.csv\n";

    std::cout << "Default configuration file created at: " << config_file_path << '\n';
    std::cout << "Please review the settings and run the simulation again." << '\n';
}

bool parse_parameters(int argc, char** argv, simulation_parameters& params) {
    using boost::is_any_of;
    using boost::split;
    using boost::token_compress_on;
    using boost::program_options::notify;
    using boost::program_options::options_description;
    using boost::program_options::parse_command_line;
    using boost::program_options::parse_config_file;
    using boost::program_options::store;
    using boost::program_options::value;
    using boost::program_options::variables_map;
    using std::ifstream;
    using std::string;
    using std::vector;

    try {
        options_description generic("Generic options");
        generic.add_options()("help,h", "Produce help message")("config,c", value<string>()->default_value("config.ini"), "Path to configuration file");

        options_description config("Configuration");
        config.add_options()("inertia_tensor", value<string>(), "[kg m^2] Spacecraft inertia tensor (9 values)")(
            "initial_angular_velocity", value<string>(), "[rad/s] Initial angular velocity (3 values)")("permanent_magnet_moment", value<string>(),
                                                                                                        "[A m^2] Permanent magnet moment (3 values)")(
            "simulation_t_start", value<double>(), "[s] Simulation start time")("simulation_t_end", value<double>(), "[s] Simulation end time")(
            "simulation_dt_initial", value<double>(), "[s] Initial time step for integrator")(
            "orbit_altitude_km", value<double>(), "[km] Circular orbit altitude")("orbit_inclination_deg", value<double>(), "[deg] Orbit inclination")(
            "output_filename", value<string>(), "Output CSV filename");

        variables_map vm;
        store(parse_command_line(argc, argv, generic), vm);
        notify(vm);

        if (vm.contains("help")) {
            std::cout << generic << "\n" << config << "\n";
            return false;
        }

        string   config_file = vm["config"].as<string>();
        ifstream ifs(config_file);

        if (!ifs) {
            if (!vm["config"].defaulted()) {
                std::cerr << "Cannot open specified config file: " << config_file << '\n';
                return false;
            }
            generate_default_config(config_file);
            return false;
        }

        store(parse_config_file(ifs, config, true), vm);
        store(parse_command_line(argc, argv, generic), vm);  // Re-parse to ensure cmdline priority
        notify(vm);

        vector<string> tokens;

        {
            auto inertia_tensor = vm["inertia_tensor"].as<string>();
            tokens.clear();
            split(tokens, inertia_tensor, is_any_of(" "), token_compress_on);
            if (tokens.size() != inertia_t::MaxSizeAtCompileTime) {
                std::cerr << "Invalid ineratia tensor" << '\n';
            } else {
                params.inertia_tensor << boost::lexical_cast<double>(tokens[0]), boost::lexical_cast<double>(tokens[1]), boost::lexical_cast<double>(tokens[2]),
                    boost::lexical_cast<double>(tokens[3]), boost::lexical_cast<double>(tokens[4]), boost::lexical_cast<double>(tokens[5]),
                    boost::lexical_cast<double>(tokens[6]), boost::lexical_cast<double>(tokens[7]), boost::lexical_cast<double>(tokens[8]);
            }
        }

        {
            auto initial_angular_velocity = vm["initial_angular_velocity"].as<string>();
            tokens.clear();
            split(tokens, initial_angular_velocity, is_any_of(" "), token_compress_on);
            if (tokens.size() != 3) {
                std::cerr << "Invalid initial angular velocity" << '\n';
            } else {
                params.initial_angular_velocity << boost::lexical_cast<double>(tokens[0]), boost::lexical_cast<double>(tokens[1]),
                    boost::lexical_cast<double>(tokens[2]);
            }
        }

        {
            auto permanent_magnet_moment = vm["permanent_magnet_moment"].as<string>();
            tokens.clear();
            split(tokens, permanent_magnet_moment, is_any_of(" "), token_compress_on);
            if (tokens.size() != 3) {
                std::cerr << "Invalid permanent magnet moment" << '\n';
            } else {
                params.permanent_magnet_moment << boost::lexical_cast<double>(tokens[0]), boost::lexical_cast<double>(tokens[1]),
                    boost::lexical_cast<double>(tokens[2]);
            }
        }

        params.t_start               = vm["simulation_t_start"].as<double>();
        params.t_end                 = vm["simulation_t_end"].as<double>();
        params.dt_initial            = vm["simulation_dt_initial"].as<double>();
        params.orbit_altitude_km     = vm["orbit_altitude_km"].as<double>();
        params.orbit_inclination_deg = vm["orbit_inclination_deg"].as<double>();
        params.output_filename       = vm["output_filename"].as<string>();
    } catch (const boost::program_options::error& e) {
        std::cerr << "Could not parse options: " << e.what() << '\n';
        return false;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return false;
    }
    return true;
}

}  // namespace aos

int main(int argc, char** argv) {
    using aos::parse_parameters;
    using aos::simulation_parameters;

    simulation_parameters params;
    if (not parse_parameters(argc, argv, params)) {
        return 1;
    }
    return 0;
}
