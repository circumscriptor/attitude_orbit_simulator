import math

def calculate_dm_dt(params, m_scalar_am, h_along_rod, dh_dt):
    """
    Calculates the magnetization derivative based on the Jiles-Atherton model.
    This function mirrors the C++ implementation for verification.
    """
    absolute_error = 1e-6
    denominator_epsilon = 1e-9

    # 1. Anhysteretic Magnetization (Man) and its derivative
    heff = h_along_rod + (params['alpha'] * m_scalar_am)
    man = 0.0
    dman_dheff = 0.0

    if abs(heff / params['a']) > absolute_error:
        x = heff / params['a']
        tanhx = math.tanh(x)
        man = params['ms'] * (1.0 / tanhx - 1.0 / x)
        dman_dheff = params['ms'] / params['a'] * (1.0 - 1.0 / (tanhx * tanhx) + 1.0 / (x * x))

    # 2. Reversible and Irreversible components of dM/dH
    sign_term = 1.0 if dh_dt >= 0.0 else -1.0
    denominator = (params['k'] * sign_term) - (params['alpha'] * (man - m_scalar_am))

    if abs(denominator) < denominator_epsilon:
        denominator = math.copysign(denominator_epsilon, denominator)

    dmirr_dh = (man - m_scalar_am) / denominator
    dmrev_dh = params['c'] * dman_dheff

    dm_dh = ((1.0 - params['c']) * dmirr_dh) + dmrev_dh

    # 3. Return dM/dt using the chain rule
    return dm_dh * dh_dt

# Parameters from the C++ Test Fixture
ja_params = {
    'ms': 1.4e5,
    'a': 2.0e3,
    'k': 1.0e3,
    'c': 0.2,
    'alpha': 1.0e-3
}

# --- Values for MagnetizationDerivativeFromH test ---
m1 = 5.0e4
h1 = 1.5e3
dh_dt1 = 1.0e2
dm_dt1 = calculate_dm_dt(ja_params, m1, h1, dh_dt1)
print(f"Result for MagnetizationDerivativeFromH: {dm_dt1:.8f}")

# --- Values for MagnetizationDerivative test ---
# First, we need to calculate h_along_rod and dh_dt from B and omega
vacuum_permeability = 4 * math.pi * 1e-7
orientation_body = (1.0, 0.0, 0.0)
b_body_t = (0.002, 0.001, 0.0)
omega_rad_s = (0.0, 0.0, 0.1)

h_along_rod2 = (b_body_t[0] * orientation_body[0]) / vacuum_permeability

# Cross product: -omega.cross(b_body_t)
cross_x = -(omega_rad_s[1] * b_body_t[2] - omega_rad_s[2] * b_body_t[1])
cross_y = -(omega_rad_s[2] * b_body_t[0] - omega_rad_s[0] * b_body_t[2])
cross_z = -(omega_rad_s[0] * b_body_t[1] - omega_rad_s[1] * b_body_t[0])
db_dt_body = (cross_x, cross_y, cross_z)
dh_dt2 = (db_dt_body[0] * orientation_body[0]) / vacuum_permeability

dm_dt2 = calculate_dm_dt(ja_params, m1, h_along_rod2, dh_dt2)
print(f"Result for MagnetizationDerivative:       {dm_dt2:.8f}")


# --- Values for NegativeDhDt test ---
m3 = 5.0e4
h3 = 1.5e3
dh_dt3 = -1.0e2
dm_dt3 = calculate_dm_dt(ja_params, m3, h3, dh_dt3)
print(f"Result for NegativeDhDt:              {dm_dt3:.8f}")
