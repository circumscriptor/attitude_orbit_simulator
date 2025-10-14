### 1. The State Vector (Y)

First, define the complete state of your system as a single vector $Y$. This vector contains all the quantities that change over time and that you need to solve for.

$$
\begin{equation}
Y = [r,v,q,\omega,M_{irr}]
\end{equation}
$$

Where:
*   $r$: Position vector of the spacecraft in an inertial frame (ECI). (3x1 vector)
*   $v$: Velocity vector of the spacecraft in an inertial frame (ECI). (3x1 vector)
*   $q$: Attitude quaternion, representing the rotation from the inertial frame to the spacecraft's body frame. (4x1 vector, $[w, x, y, z]$)
*   $\omega$: Angular velocity vector of the spacecraft in its own body frame. (3x1 vector)
*   $M_{irr}$: A vector containing the scalar irreversible magnetization for *each* of the $N$ hysteresis rods. ($[M_irr,1, M_irr,2, ..., M_irr,N]$)

### 2. The System of First-Order Differential Equations $\left(\frac{dY}{dt}\right)$

Function that takes the current time $t$ and the current state vector $Y$ as input and returns the time derivative of every component of $Y$.

$$
\begin{equation}
\frac{d\mathbf{Y}}{dt} =
\begin{bmatrix}
\dot{\mathbf{r}} \\
\dot{\mathbf{v}} \\
\dot{\mathbf{q}} \\
\dot{\boldsymbol{\omega}} \\
\dot{\mathbf{M}}_{irr}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{v} \\
-\frac{G M_{Earth}}{||\mathbf{r}||^3} \mathbf{r} \\
\frac{1}{2} \mathbf{q} \otimes [0, \omega_x, \omega_y, \omega_z] \\
\mathbf{I}^{-1} \left( \boldsymbol{\tau}_{total} - \boldsymbol{\omega} \times (\mathbf{I}\boldsymbol{\omega}) \right) \\
\left[ \frac{dM_{irr,1}}{dt}, \frac{dM_{irr,2}}{dt}, \dots, \frac{dM_{irr,N}}{dt} \right]^T
\end{bmatrix}
\end{equation}
$$

### 3. Calculation Workflow (The $f(t, Y)$ function)

Perform a sequence of calculations using the current state $Y$ and system parameters.

**Step 1: Unpack the State Vector $Y$**
   *   Extract $r$, $v$, $q$, $\omega$, and the array of $M_{irr,i}$ values.
   *   Normalize the quaternion $q$ to ensure it remains a unit quaternion: $q = \frac{q}{||q||}$.

**Step 2: Calculate Intermediate Environmental Variables**
*   **$B_{body}$ (Magnetic Field in Body Frame):**
    1.  Use $r$ and $t$ to call a geomagnetic model (e.g., WMM) to get the B-field in an Earth-fixed frame ($B_{ECEF}$).
    2.  Perform the coordinate transformation from the Earth-fixed frame to the inertial frame ($B_{ECI}$).
    3.  Use the attitude quaternion $q$ to rotate $B_{ECI}$ into the spacecraft's body frame to get $B_{body}$.
* **$\frac{dH_i}{dt}$ (Rate of change of field along each rod $i$):**
    1.  Calculate the rotational component:
        $$
        \begin{equation}
        \left(\frac{dB}{dt}\right)_{rot} = \omega \times B_{body}
        \end{equation}
        $$
    2.  Calculate the orbital component $(\frac{dB}{dt})_{orb}$ by numerically differencing the B-field along the orbit (e.g., using the current velocity $v$ to find the position at $t+\delta{t}$).
    3.  Combine them:
        $$
        \begin{equation}
        \frac{dH_i}{dt} =
        \left(\frac{1}{\mu_0}\right) * \left[ \left( \left(\frac{dB}{dt}\right)_{orb} + \left(\frac{dB}{dt}\right)_{rot} \right) ⋅ b_i \right]
        \end{equation}
        $$

**Step 3: Calculate the Hysteresis State and Total Magnetization $M_{total,i}$ (for each rod $i$)**

$$
\begin{align}
    H_i &= (B_{body} ⋅ b_i) / \mu_0 \\
    H_{eff,i} &= H_i + α * M_{total,i_{prev}} \\
    M_{an,i} &= M_s * [\coth(H_{eff,i} / a) - a / H_{eff,i}] \\
    \frac{dM_{irr,i}}{dt} &= \frac{M_{an,i} - M_{irr,i}}{k\delta} \cdot \frac{dH_i}{dt} \\
    M_{total,i} &= M_irr,i + c * (M_{an,i} - M_irr,i)
\end{align}
$$

**Step 4: Calculate Torques**

$$
\begin{align}
\tau_p &= m_p \times B_{body} \\
\mu_{h,i} &= M_{total,i} \cdot V_i \cdot b_i \\
\tau_{h,i} &= \mu_{h,i} \times B_{body} \\
\tau_{total} &= \tau_p + \sum_{i=0}^{N}{\tau_{h,i}}
\end{align}
$$

**Step 5: Assemble the Final Derivative Vector $dY/dt$**
$$
\begin{align}
    \frac{dr}{dt} &= v \\
    \frac{dv}{dt} &= -G \cdot M_{earth} \cdot \frac{r}{||r||^3} \\
    \frac{dq}{dt} &= 0.5 \cdot q \otimes [0, \omega] \\
    \frac{d\omega}{dt} &= I^{-1} \cdot (\tau_{total} - \omega \times (I\omega)) \\
    \frac{dM_{irr}}{dt} &= \text{The vector of } \frac{dM_{irr,i}}{dt} \text{ calculated in Step 3.}
\end{align}
$$
