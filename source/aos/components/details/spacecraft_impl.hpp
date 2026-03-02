#pragma once

#include "aos/components/hysteresis_rods.hpp"
#include "aos/components/inertia_tensor.hpp"
#include "aos/components/permanent_magnet.hpp"
#include "aos/components/spacecraft.hpp"
#include "aos/components/spacecraft_faces.hpp"
#include "aos/core/state.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

namespace aos {

class spacecraft_impl : public spacecraft {
public:

    spacecraft_impl(const spacecraft_impl&)                    = delete;
    spacecraft_impl(spacecraft_impl&&)                         = delete;
    auto operator=(const spacecraft_impl&) -> spacecraft_impl& = delete;
    auto operator=(spacecraft_impl&&) -> spacecraft_impl&      = delete;

    explicit spacecraft_impl(const spacecraft_properties& properties);
    ~spacecraft_impl() override;

    [[nodiscard]] auto mass_kg() const -> double override;
    [[nodiscard]] auto inertia() const -> const inertia_tensor& override;
    [[nodiscard]] auto faces() const -> const spacecraft_faces& override;
    [[nodiscard]] auto magnet() const -> const permanent_magnet& override;
    [[nodiscard]] auto hystresis() const -> const hysteresis_rods& override;

    void derivative(const environment_effects& env, const system_state& current_state, system_state& state_derivative) const override;

    // sums permanent magnet, gyroscopic, and gravity gradient torques
    [[nodiscard]] auto compute_torques(const vec3& omega, const vec3& b_body, const vec3& r_body, double earth_mu) const -> vec3;

private:

    double           _mass_kg;  // [kg] Mass
    inertia_tensor   _inertia;
    spacecraft_faces _faces;
    permanent_magnet _magnet;
    hysteresis_rods  _hystresis;
};

}  // namespace aos
