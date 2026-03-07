#pragma once

#include "aos/components/hysteresis_rod.hpp"
#include "aos/components/hysteresis_rods.hpp"
#include "aos/components/inertia_tensor.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/components/spacecraft_faces.hpp"
#include "aos/components/spacecraft_shape.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

#include <memory>

namespace aos {

struct spacecraft_properties {
    real                        mass_kg;
    spacecraft_shape            shape;
    hysteresis_parameters       hysteresis;
    permanent_magnet_properties magnet;
    hysteresis_rods_properties  rods;

    void from_toml(const toml_table& table);
    void debug_print() const;
};

class spacecraft {
public:

    spacecraft(const spacecraft&)                    = delete;
    spacecraft(spacecraft&&)                         = delete;
    auto operator=(const spacecraft&) -> spacecraft& = delete;
    auto operator=(spacecraft&&) -> spacecraft&      = delete;

    explicit spacecraft(const spacecraft_properties& properties);
    ~spacecraft();

    [[nodiscard]] auto mass_kg() const -> double;
    [[nodiscard]] auto inertia() const -> const inertia_tensor&;
    [[nodiscard]] auto faces() const -> const spacecraft_faces&;
    [[nodiscard]] auto magnet() const -> const permanent_magnet&;
    [[nodiscard]] auto hystresis() const -> const hysteresis_rods&;

    void derivative(const environment_effects& env, const system_state& current_state, system_state& state_derivative) const;

    static auto create(const spacecraft_properties& properties) -> std::shared_ptr<spacecraft>;

protected:

    // sums permanent magnet, gyroscopic, and gravity gradient torques
    [[nodiscard]] auto compute_torques(const vec3& omega, const vec3& b_body, const vec3& r_body, real earth_mu) const -> vec3;

private:

    real             _mass_kg;  // [kg] Mass
    inertia_tensor   _inertia;
    spacecraft_faces _faces;
    permanent_magnet _magnet;
    hysteresis_rods  _hystresis;
};

}  // namespace aos
