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
    double                      mass_kg;
    spacecraft_shape            shape;
    hysteresis_parameters       hysteresis;
    permanent_magnet_properties magnet;
    hysteresis_rods_properties  rods;

    void from_toml(const toml_table& table);
    void debug_print() const;
};

class spacecraft {
public:

    spacecraft()                                     = default;
    spacecraft(const spacecraft&)                    = delete;
    spacecraft(spacecraft&&)                         = delete;
    auto operator=(const spacecraft&) -> spacecraft& = delete;
    auto operator=(spacecraft&&) -> spacecraft&      = delete;

    virtual ~spacecraft();

    [[nodiscard]] virtual auto mass_kg() const -> double                   = 0;
    [[nodiscard]] virtual auto inertia() const -> const inertia_tensor&    = 0;
    [[nodiscard]] virtual auto faces() const -> const spacecraft_faces&    = 0;
    [[nodiscard]] virtual auto magnet() const -> const permanent_magnet&   = 0;
    [[nodiscard]] virtual auto hystresis() const -> const hysteresis_rods& = 0;

    virtual void derivative(const environment_effects& env, const system_state& current_state, system_state& state_derivative) const = 0;

    static auto create(const spacecraft_properties& properties) -> std::shared_ptr<spacecraft>;
};

}  // namespace aos
