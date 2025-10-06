### Preamble: Definition of Variables

The state of the system at any time `t` is defined by the following variables:
*   **Attitude Quaternion:** $\mathbf{q}(t) = [q_w, q_x, q_y, q_z]^T$
*   **Angular Velocity:** $\boldsymbol{\omega}(t) = [\omega_x, \omega_y, \omega_z]^T$
*   **Hysteresis Rod Magnetizations:** $M_i(t)$ for each of the $N$ rods.

The evolution of these states depends on the following constant parameters and time-varying inputs.

**Constant Parameters:**
*   **Inertia Tensor (diagonal):** $I_{xx}, I_{yy}, I_{zz}$
*   **Permanent Magnet Moment:** $[m_{px}, m_{py}, m_{pz}]^T$
*   **Hysteresis Rod Orientations:** $[b_{ix}, b_{iy}, b_{iz}]^T$ for each rod $i$.
*   **Hysteresis Rod Volume:** $V_{h,i}$ for each rod $i$.
*   **Jiles-Atherton Parameters:** $M_{s,i}, a_i, k_i, c_i, \alpha_i$ for each rod $i$.
*   **Vacuum Permeability:** $\mu_0$

**Time-Varying Inputs:**
*   **Inertial Magnetic Field:** $\mathbf{B}_{in}(t) = [B_{ix}(t), B_{iy}(t), B_{iz}(t)]^T$, provided by the environment model.

---

### I. The 4 Attitude Kinematics Equations

These equations describe how the attitude quaternion evolves based on the angular velocity.

$$
\frac{dq_w}{dt} = \frac{1}{2} (- \omega_x q_x - \omega_y q_y - \omega_z q_z)
$$

$$
\frac{dq_x}{dt} = \frac{1}{2} ( \omega_x q_w + \omega_z q_y - \omega_y q_z)
$$

$$
\frac{dq_y}{dt} = \frac{1}{2} ( \omega_y q_w - \omega_z q_x + \omega_x q_z)
$$

$$
\frac{dq_z}{dt} = \frac{1}{2} ( \omega_z q_w + \omega_y q_x - \omega_x q_y)
$$

---

### II. The 3 Rotational Dynamics Equations

These are the unrolled Euler's equations for angular acceleration. To write them fully, we must first define the components of the magnetic field in the body frame, $\mathbf{B}_{body} = [B_x, B_y, B_z]^T$.

**A. Transformation of the Magnetic Field**

The magnetic field is transformed from the inertial frame to the body frame using the transpose of the quaternion's rotation matrix:

$$
\begin{bmatrix} B_x \\ B_y \\ B_z \end{bmatrix}
=
\begin{bmatrix}
1 - 2(q_y^2 + q_z^2) & 2(q_x q_y + q_w q_z) & 2(q_x q_z - q_w q_y) \\
2(q_x q_y - q_w q_z) & 1 - 2(q_x^2 + q_z^2) & 2(q_y q_z + q_w q_x) \\
2(q_x q_z + q_w q_y) & 2(q_y q_z - q_w q_x) & 1 - 2(q_x^2 + q_y^2)
\end{bmatrix}
\begin{bmatrix} B_{ix} \\ B_{iy} \\ B_{iz} \end{bmatrix}
$$

**B. The Dynamics Equations**

Using the $B_x, B_y, B_z$ calculated above:

$$
\frac{d\omega_x}{dt} = \frac{1}{I_{xx}} \left[
(m_{py} B_z - m_{pz} B_y)
+ \sum_{i=1}^{N} \frac{V_{h,i} M_i(t)}{\mu_0} (b_{iy} B_z - b_{iz} B_y)
+ (I_{yy} - I_{zz}) \omega_y \omega_z
\right]
$$

$$
\frac{d\omega_y}{dt} = \frac{1}{I_{yy}} \left[
(m_{pz} B_x - m_{px} B_z)
+ \sum_{i=1}^{N} \frac{V_{h,i} M_i(t)}{\mu_0} (b_{iz} B_x - b_{ix} B_z)
+ (I_{zz} - I_{xx}) \omega_z \omega_x
\right]
$$

$$
\frac{d\omega_z}{dt} = \frac{1}{I_{zz}} \left[
(m_{px} B_y - m_{py} B_x)
+ \sum_{i=1}^{N} \frac{V_{h,i} M_i(t)}{\mu_0} (b_{ix} B_y - b_{iy} B_x)
+ (I_{xx} - I_{yy}) \omega_x \omega_y
\right]
$$

---

### III. The N Hysteresis Rod Dynamics Equations

For each of the $N$ hysteresis rods (indexed by $i$), there is a differential equation for its scalar magnetization $M_i$. This is the unrolled Jiles-Atherton model.

**A. Intermediate Variables for Rod $i$**

First, we define the magnetic field along the rod's axis ($H_i$) and its rate of change ($dH_i/dt$):

$$
H_i = \frac{1}{\mu_0} (B_x b_{ix} + B_y b_{iy} + B_z b_{iz})
$$

$$
\frac{dH_i}{dt} = \frac{1}{\mu_0} \left[
  ( \omega_z B_y - \omega_y B_z ) b_{ix} +
  ( \omega_x B_z - \omega_z B_x ) b_{iy} +
  ( \omega_y B_x - \omega_x B_y ) b_{iz}
\right]
$$

**B. The Jiles-Atherton Differential Equation for Rod $i$**

The rate of change of magnetization for rod $i$ is given by:

$$
\frac{dM_i}{dt} = \left[
(1 - c_i) \frac{M_{an,i} - M_i}{k_i \cdot \operatorname{sign}(\frac{dH_i}{dt})}
+
c_i \left( \frac{M_{s,i}}{a_i} \left(1 - \coth^2\left(\frac{H_{e,i}}{a_i}\right) + \frac{a_i^2}{H_{e,i}^2}\right) \right)
\right] \times \frac{dH_i}{dt}
$$

Where the terms $M_{an,i}$ and $H_{e,i}$ are themselves functions of the current state:

$$
H_{e,i} = H_i + \alpha_i M_i
$$

$$
M_{an,i} = M_{s,i} \left( \coth\left(\frac{H_{e,i}}{a_i}\right) - \frac{a_i}{H_{e,i}} \right)
$$
