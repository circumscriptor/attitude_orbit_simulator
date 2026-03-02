#pragma once

#include "aos/components/spacecraft_face.hpp"
#include "aos/components/spacecraft_shape.hpp"
#include "aos/core/types.hpp"
#include "aos/environment/environment.hpp"

namespace aos {

class spacecraft_faces {
public:

    spacecraft_faces(const spacecraft_faces&)                    = delete;
    spacecraft_faces(spacecraft_faces&&)                         = delete;
    auto operator=(const spacecraft_faces&) -> spacecraft_faces& = delete;
    auto operator=(spacecraft_faces&&) -> spacecraft_faces&      = delete;

    explicit spacecraft_faces(const spacecraft_shape& shape);
    explicit spacecraft_faces(const spacecraft_custom& shape);
    explicit spacecraft_faces(const spacecraft_uniform& shape);
    ~spacecraft_faces() = default;

    // compute drag and srp torque+force
    [[nodiscard]] auto compute_face_effects(const environment_effects& env, const quat& q_att, const quat& q_inv, const vec3& omega_body) const -> face_effects;

    // compute drag and srp torque+force
    [[nodiscard]] auto compute_faces_effects_with_forces(const environment_effects& env, const quat& q_inv, const vec3& omega_body) const
        -> face_effects_with_forces;

protected:

    void uniform(const spacecraft_uniform& shape);

private:

    spacecraft_face_array _faces;
};

}  // namespace aos
