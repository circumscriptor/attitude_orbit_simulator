#include "aos/core/types.hpp"
#include "aos/simulation/config.hpp"
#include "aos/simulation/simulation.hpp"
#include "aos/verify/hysteresis.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/constants.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/detail/parsers.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void parse_option_vec3(std::vector<std::string>& tokens, const std::string& input, aos::vec3& output) {
    tokens.clear();
    tokens.reserve(3);

    boost::split(tokens, input, boost::is_any_of(","), boost::token_compress_on);
    if (tokens.size() != 3) {
        throw std::runtime_error("Invalid option size, expceted 3");
    }

    output <<                                    //
        boost::lexical_cast<double>(tokens[0]),  //
        boost::lexical_cast<double>(tokens[1]),  //
        boost::lexical_cast<double>(tokens[2]);  //
}

}  // namespace

int main(int argc, char** argv) {
    using boost::program_options::bool_switch;
    using boost::program_options::error;
    using boost::program_options::notify;
    using boost::program_options::options_description;
    using boost::program_options::parse_command_line;
    using boost::program_options::store;
    using boost::program_options::value;
    using boost::program_options::variables_map;

    static const unsigned help_line_width = 40;

    auto params = aos::simulation_parameters::get_default();

    options_description generic("Generic options");
    generic.add_options()                                                                //
        ("help,h", "Print this help message")                                            //
        ("output,o", value<std::string>()->default_value("output.csv"), "Output file");  //

    options_description simulation_parameters("Simulation parameters");
    simulation_parameters.add_options()                                                                              //
        ("mass", value<double>(&params.spacecraft.mass_g), "Spacecraft mass [g]")                                    //
        ("width", value<double>(), "Spacecraft width [m]")                                                           //
        ("height", value<double>(), "Spacecraft height [m]")                                                         //
        ("length", value<double>(), "Spacecraft length [m]")                                                         //
        ("magnet-remanence", value<double>(&params.spacecraft.magnet_remanence), "Permanent magnet remanence")       //
        ("magnet-length", value<double>(&params.spacecraft.magnet_length), "Permanent magnet length [m]")            //
        ("magnet-diameter", value<double>(&params.spacecraft.magnet_diameter), "Permanent magnet diameter [m]")      //
        ("no-rods", "Do not use hysteresis rods")                                                                    //
        ("rod-volume", value<double>(&params.spacecraft.hysteresis_rod_volume), "Volume of hysteresis rod in [m3]")  //
        ("rod-orientation", value<std::vector<std::string>>(), "Hysteresis rod orientation (multiple)")              //
        ("altitude", value<double>(&params.environment.orbit_altitude_km), "Orbit altitude [km]")                    //
        ("inclination", value<double>(&params.environment.orbit_inclination_deg), "Orbit inclination [deg]")         //
        ("angular-velocity", value<std::string>(), "Initial angular velocity xyz [rad/s]")                           //
        ("t-start", value<double>(&params.t_start), "Simulation start time [s]")                                     //
        ("t-end", value<double>(&params.t_end), "Simulation end time [s]")                                           //
        ("dt", value<double>(&params.dt_initial), "Initial simulation time step [s]");                               //

    options_description hysteresis_parameters("Hysteresis parameters");
    hysteresis_parameters.add_options()                                                                                //
        ("hysteresis-ms", value<double>(&params.spacecraft.hysteresis_params.ms), "Saturation Magnetization [A/m]")    //
        ("hysteresis-a", value<double>(&params.spacecraft.hysteresis_params.a), "Anhysteretic shape parameter [A/m]")  //
        ("hysteresis-k", value<double>(&params.spacecraft.hysteresis_params.k), "Pinning energy density [A/m]")        //
        ("hysteresis-c", value<double>(&params.spacecraft.hysteresis_params.c), "Reversibility coefficient [0-1]")     //
        ("hysteresis-alpha", value<double>(&params.spacecraft.hysteresis_params.alpha), "Inter-domain coupling");      //

    options_description other("Other options");
    other.add_options()                                                                                                                        //
        ("higher-order", bool_switch(&params.higher_order), "Use higher order solver")                                                         //
        ("absolute-error", value<double>(&params.absolute_error), "Integration solver's absolute error")                                       //
        ("relative-error", value<double>(&params.relative_error), "Integration solver's relative error")                                       //
        ("verify-hysteresis", "Calculate hysteresis curve for the given material instead of simulation")                                       //
        ("no-observe-element", bool_switch(&params.observer.exclude_elements), "Exclude per-element values from output")                       //
        ("no-observe-magnitude", bool_switch(&params.observer.exclude_magnitudes), "Exclude magnitude values from output")                     //
        ("hysteresis-material", value<std::string>()->default_value("hymu80"), "The material for hysteresis generation (not supported yet)");  //

    options_description options(help_line_width);
    options.add(generic).add(simulation_parameters).add(hysteresis_parameters).add(other);

    variables_map vm;
    try {
        store(parse_command_line(argc, argv, options), vm);
        notify(vm);
    } catch (const error& err) {
        std::cerr << "Could not parse command line options: " << err.what() << '\n';
        return 1;
    } catch (const std::exception& ex) {
        std::cerr << "Could not parse command line options: " << ex.what() << '\n';
        return 1;
    }

    if (vm.contains("help")) {
        std::cout << "Attitude Orbit Simulator for Passive AOCS\n";
        std::cout << options << '\n';
        return 0;
    }

    // Parse remaining params here
    try {
        std::vector<std::string> tokens;

        if (vm.contains("width")) {
            params.spacecraft.dim_m.x() = vm["width"].as<double>();
        }

        if (vm.contains("height")) {
            params.spacecraft.dim_m.y() = vm["height"].as<double>();
        }

        if (vm.contains("length")) {
            params.spacecraft.dim_m.z() = vm["length"].as<double>();
        }

        if (vm.contains("no-rods")) {
            params.spacecraft.hysteresis_rod_orientations.clear();
        } else if (vm.contains("rod-orientation")) {
            params.spacecraft.hysteresis_rod_orientations.clear();

            const auto& orientations = vm["rod-orientation"].as<std::vector<std::string>>();
            for (const auto& orientation : orientations) {
                parse_option_vec3(tokens, orientation, params.spacecraft.hysteresis_rod_orientations.emplace_back());
            }
        }

        if (vm.contains("angular-velocity")) {
            const auto& angular_velocity = vm["angular-velocity"].as<std::string>();
            parse_option_vec3(tokens, angular_velocity, params.angular_velocity);
        }
    } catch (const error& err) {
        std::cerr << "Could not parse parameters: " << err.what() << '\n';
        return 1;
    } catch (const std::exception& ex) {
        std::cerr << "Could not parse parameters: " << ex.what() << '\n';
        return 1;
    }

    try {
        if (vm.contains("verify-hysteresis")) {
            params.spacecraft.hysteresis_params.debug_print();
            aos::verify_hysteresis(vm["output"].as<std::string>(), params.spacecraft.hysteresis_params);
        } else {
            params.debug_print();
            aos::run_simulation(vm["output"].as<std::string>(), params);
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
    return 0;
}
