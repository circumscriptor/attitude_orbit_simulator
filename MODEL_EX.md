### 1. The State Vector $\mathbf{Y}$

$$
\mathbf{Y} = [\mathbf{r}, \mathbf{v}, \mathbf{q}, \boldsymbol{\omega}, \mathbf{M}_{irr}]^T
$$

*   $\mathbf{r}$: Position in **ECI** ($m$).
*   $\mathbf{v}$: Velocity in **ECI** ($m/s$).
*   $\mathbf{q}$: Quaternion (ECI $\to$ Body) ($[w, x, y, z]$, unit norm).
*   $\boldsymbol{\omega}$: Angular velocity in **Body** ($rad/s$).
*   $\mathbf{M}_{irr}$: Irreversible magnetization of hysteresis rods ($A/m$).

---

### 2. The Differential Equation $\frac{d\mathbf{Y}}{dt}$

$$
\frac{d\mathbf{Y}}{dt} = \begin{bmatrix}
\mathbf{v} \\
-\frac{\mu}{\lVert\mathbf{r}\rVert^3}\mathbf{r} + \mathbf{a}_{pert, ECI} \\
\frac{1}{2} \mathbf{q} \otimes [0, \boldsymbol{\omega}]^T \\
\mathbf{I}^{-1} \left( \boldsymbol{\tau}_{mag} + \boldsymbol{\tau}_{aero} + \boldsymbol{\tau}_{grav} - \boldsymbol{\omega} \times (\mathbf{I}\boldsymbol{\omega}) \right) \\
\dot{\mathbf{M}}_{irr}
\end{bmatrix}
$$

*Note: $\mu$ is retrieved via `GeographicLib::GravityModel::GM()`.*

---

### 3. Simulation Workflow $f(t, \mathbf{Y})$

#### Step 1: Time and Frame Setup
We must bridge the gap between the Inertial frame (ECI) and the Earth-Fixed frame (ECEF/Geodetic).

1.  **Calculate GMST (Greenwich Mean Sidereal Time):**
    Compute angle $\theta_{gst}$ based on Julian Date of time $t$.
2.  **Construct Rotation Matrix $\mathbf{R}_{ECEF}^{ECI}$:**
    $$
    \mathbf{R}_{ECEF}^{ECI} = \begin{bmatrix}
    \cos\theta & -\sin\theta & 0 \\
    \sin\theta & \cos\theta & 0 \\
    0 & 0 & 1
    \end{bmatrix}
    $$
3.  **Position Conversion:**
    $$ \mathbf{r}_{ECEF} = (\mathbf{R}_{ECEF}^{ECI})^T \mathbf{r}_{ECI} $$

#### Step 2: GeographicLib Coordinates
Use `Geocentric` to get Geodetic coordinates and the local frame rotation.

```cpp
GeographicLib::Geocentric earth(Constants::WGS84_a(), Constants::WGS84_f());
double lat, lon, h;
// Computes lat/lon/h AND the rotation from ENU to ECEF (M)
earth.Reverse(r_ecef.x, r_ecef.y, r_ecef.z, lat, lon, h, M);
```
*   **Input:** $\mathbf{r}_{ECEF}$.
*   **Output:** Latitude ($\phi$), Longitude ($\lambda$), Height ($h$).
*   **Output:** $\mathbf{R}_{ENU}^{ECEF}$ (Returned by `Reverse` as a $3\times3$ matrix or computed via `LocalCartesian`).

#### Step 3: Environmental Vectors

**A. Gravity Perturbations (GravityModel)**
1.  Get disturbance in Local Tangent Plane (ENU).
    ```cpp
    grav.Disturbance(lat, lon, h, gx, gy, gz); // gx=East, gy=North, gz=Up
    Vector3d a_pert_ENU(gx, gy, gz);
    ```
2.  Transform to Inertial Frame (ECI) for the integrator:
    $$ \mathbf{a}_{pert, ECI} = \mathbf{R}_{ECEF}^{ECI} \cdot \mathbf{R}_{ENU}^{ECEF} \cdot \mathbf{a}_{pert, ENU} $$

**B. Magnetic Field (MagneticModel)**
1.  Get field in ENU.
    ```cpp
    mag(year, lat, lon, h, Bx, By, Bz); // nT
    Vector3d B_ENU(Bx, By, Bz);
    B_ENU *= 1e-9; // Convert nT to Tesla
    ```
2.  Transform to Body Frame (for Torque & Hysteresis):
    *   First to ECI: $\mathbf{B}_{ECI} = \mathbf{R}_{ECEF}^{ECI} \cdot \mathbf{R}_{ENU}^{ECEF} \cdot \mathbf{B}_{ENU}$
    *   Then to Body: $\mathbf{B}_{body} = \mathbf{R}_{ECI}^{Body}(\mathbf{q}) \cdot \mathbf{B}_{ECI}$

#### Step 4: Magnetic Field Derivative ($\dot{\mathbf{B}}$)
Hysteresis depends on how fast the field changes *inside the rod*. This has two components: Orbit movement and Satellite tumbling.

$$
\frac{d\mathbf{B}_{body}}{dt} = \underbrace{\mathbf{R}_{ECI}^{Body}(\mathbf{q}) \frac{\mathbf{B}_{ECI}(t) - \mathbf{B}_{ECI}(t-\Delta t)}{\Delta t}}_{\text{Orbital Change}} - \underbrace{\boldsymbol{\omega} \times \mathbf{B}_{body}}_{\text{Rotational Change}}
$$

*   **Logic:**
    *   If $t=0$, set $\dot{\mathbf{B}}_{body} = 0$ (or assume a small initial rate).
    *   Store $\mathbf{B}_{ECI}(t)$ to use as $\mathbf{B}_{ECI}(t-\Delta t)$ in the *next* step.
    *   The $-\boldsymbol{\omega} \times \mathbf{B}$ term is the transport theorem (derivative of a vector in a rotating frame).

#### Step 5: Hysteresis Dynamics ($\dot{\mathbf{M}}_{irr}$)
For each rod $i$ aligned with body axis $\mathbf{u}_i$:

1.  **Project Field and Rate onto Rod:**
    $$
    H_i = \frac{1}{\mu_0} (\mathbf{B}_{body} \cdot \mathbf{u}_i)
    $$
    $$
    \dot{H}_i = \frac{1}{\mu_0} \left( \frac{d\mathbf{B}_{body}}{dt} \cdot \mathbf{u}_i \right)
    $$

2.  **Calculate Theoretical Anysteretic Magnetization ($M_{an}$):**
    $$ H_{eff} = H_i + \alpha M_{irr, i} $$
    $$ M_{an} = M_s \left( \coth\left(\frac{H_{eff}}{a}\right) - \frac{a}{H_{eff}} \right) $$

3.  **Calculate Derivative (The Hysteresis Differential Equation):**
    Note: The classic equation gives $\frac{dM}{dH}$. We need $\frac{dM}{dt}$.
    $$
    \frac{dM_{irr, i}}{dH} = \frac{M_{an} - M_{irr, i}}{k \delta} \quad \text{where } \delta = \text{sign}(\dot{H}_i)
    $$
    **Chain Rule:**
    $$
    \frac{dM_{irr, i}}{dt} = \left( \frac{M_{an} - M_{irr, i}}{k \cdot \text{sign}(\dot{H}_i)} \right) \cdot \dot{H}_i
    $$

#### Step 6: Torque Summation
Calculate the total dipole moment $\mathbf{m}_{total}$ (Permanent + Hysteresis rods).

$$
\mathbf{m}_{rods} = \sum_{i=1}^{N} (V_{rod} \cdot (M_{irr, i} + \chi_r H_i) \cdot \mathbf{u}_i)
$$
*(Note: $\chi_r H_i$ adds the reversible/linear component if explicit in your model, otherwise just $M_{total}$)*

$$
\boldsymbol{\tau}_{mag} = (\mathbf{m}_{perm} + \mathbf{m}_{rods}) \times \mathbf{B}_{body}
$$

Use this $\boldsymbol{\tau}_{mag}$ in the final derivative vector defined in Section 2.
