### 1. The State Vector $\mathbf{Y}$

First, define the complete state of your system as a single vector $Y$. This vector contains all the quantities that change over time and that you need to solve for.

$$
\begin{equation}
\mathbf{Y} = [\mathbf{r},\mathbf{v},\mathbf{q},\boldsymbol{\omega},M_{irr}]
\end{equation}
$$

Where:
*   $\mathbf{r}$ - Position vector of the spacecraft in an inertial frame (ECI). ($3 \times 1$ vector)
*   $\mathbf{v}$ - Velocity vector of the spacecraft in an inertial frame (ECI). ($3 \times 1$ vector)
*   $\mathbf{q}$ - Attitude quaternion, representing the rotation from the inertial frame to the spacecraft's body frame. ($4 \times 1$ vector, $[w, x, y, z]$)
*   $\boldsymbol{\omega}$ - Angular velocity vector of the spacecraft in its own body frame. ($3 \times 1$ vector)
*   $M_{irr}$ - A vector containing the scalar irreversible magnetization for *each* of the $N$ hysteresis rods. ($[M_irr,1, M_irr,2, ..., M_irr,N]$)

### 2. The System of First-Order Differential Equations $\frac{d\mathbf{Y}}{dt}$

Function that takes the current time $t$ and the current state vector $Y$ as input and returns the time derivative of every component of $Y$.

$$
\begin{equation}
\frac{d\mathbf{Y}}{dt} = \begin{bmatrix}
\dot{\mathbf{r}} \\
\dot{\mathbf{v}} \\
\dot{\mathbf{q}} \\
\dot{\boldsymbol{\omega}} \\
\dot{\mathbf{M}}_{irr}
\end{bmatrix} = \begin{bmatrix}
\mathbf{v} \\
-\frac{G M_{Earth}}{\lVert\mathbf{r}\rVert^3} \mathbf{r} \\
\frac{1}{2} \mathbf{q} \otimes [0, \omega_x, \omega_y, \omega_z] \\
\mathbf{I}^{-1} \left( \boldsymbol{\tau}_{total} - \boldsymbol{\omega} \times (\mathbf{I}\boldsymbol{\omega}) \right) \\
\left[ \frac{dM_{irr,1}}{dt}, \frac{dM_{irr,2}}{dt}, \dots, \frac{dM_{irr,N}}{dt} \right]^T
\end{bmatrix}
\end{equation}
$$

### 3. Calculation Workflow $f(t, \mathbf{Y})$

Perform a sequence of calculations using the current state $Y$ and system parameters.

**Step 1: Unpack the State Vector $Y$**
   *   Extract $\mathbf{r}$, $\mathbf{v}$, $\mathbf{q}$, $\boldsymbol{\omega}$, and the array of $M_{irr,i}$ values.
   *   Normalize the quaternion $\mathbf{q}$ to ensure it remains a unit quaternion:

$$
\begin{equation}
\mathbf{q} = \frac{\mathbf{q}}{\lVert\mathbf{q}\rVert}
\end{equation}
$$

**Step 2: Calculate Intermediate Environmental Variables**
*   **$B_{body}$ (Magnetic Field in Body Frame):**
    1.  Use $\mathbf{r}$ and $t$ to call a geomagnetic model (e.g., WMM) to get the B-field in an Earth-fixed frame ($B_{ECEF}$).
    2.  Perform the coordinate transformation from the Earth-fixed frame to the inertial frame ($B_{ECI}$).
    3.  Use the attitude quaternion $\mathbf{q}$ to rotate $B_{ECI}$ into the spacecraft's body frame to get $B_{body}$.
* **$\frac{dH_i}{dt}$ (Rate of change of field along each rod $i$):**

$$
\begin{align}
\left(\frac{dB}{dt}\right)_{rot,body} &= \omega \times B_{body} \\
\left(\frac{d\mathbf{B}}{dt}\right)_{orb,eci} &= (\mathbf{v}_{eci} \cdot \nabla)\mathbf{B}_{eci} \approx \frac{\mathbf{B}_{eci}(\mathbf{r}(t) + \mathbf{v}(t)\delta t) - \mathbf{B}_{eci}(\mathbf{r}(t))}{\delta t} \\
\left(\frac{d\mathbf{B}}{dt}\right)_{orb,body} &= C(\mathbf{q}) \left(\frac{d\mathbf{B}}{dt}\right)_{orb,eci} \\
\frac{dH_i}{dt} &= \frac{1}{\mu_0} \left[ \left( \left(\frac{d\mathbf{B}}{dt}\right)_{orb, body} + \left(\frac{d\mathbf{B}}{dt}\right)_{rot,body} \right) \cdot \mathbf{b}_i \right]
\end{align}
$$

$$
$$

$$
\begin{equation}
\end{equation}
$$

**Step 3: Calculate the Hysteresis State and Total Magnetization $M_{total,i}$ (for each rod $i$)**

$$
\begin{align}
    H_i &= \frac{\mathbf{B}_{body} \cdot \mathbf{b}_i}{\mu_0} \\
    H_{eff,i} &= H_i + \alpha M_{total,i_{prev}} \\
    M_{an,i} &= M_s \left[ \coth\left(\frac{H_{eff,i}}{a}\right) - \frac{a}{H_{eff,i}} \right] \\
    \frac{dM_{irr,i}}{dt} &= \frac{M_{an,i} - M_{irr,i}}{k\delta} \frac{dH_i}{dt} \\
    M_{total,i} &= M_{irr,i} + c (M_{an,i} - M_{irr,i})
\end{align}
$$

**Step 4: Calculate Torques**

$$
\begin{align}
\boldsymbol{\tau}_p &= \mathbf{m}_p \times \mathbf{B}_{body} \\
\boldsymbol{\mu}_{h,i} &= M_{total,i} V_i \mathbf{b}_i \\
\boldsymbol{\tau}_{h,i} &= \boldsymbol{\mu}_{h,i} \times \mathbf{B}_{body} \\
\boldsymbol{\tau}_{total} &= \boldsymbol{\tau}_p + \sum_{i=1}^{N}{\boldsymbol{\tau}_{h,i}}
\end{align}
$$

**Step 5: Assemble the Final Derivative Vector $\frac{dY}{dt}$**

$$
\begin{align}
    \frac{d\mathbf{r}}{dt} &= \mathbf{v} \\
    \frac{d\mathbf{v}}{dt} &= -\frac{G M_{earth}}{\lVert\mathbf{r}\rVert^3} \mathbf{r} \\
    \frac{d\mathbf{q}}{dt} &= \frac{1}{2} \mathbf{q} \otimes [0, \boldsymbol{\omega}] \\
    \frac{d\boldsymbol{\omega}}{dt} &= \mathbf{I}^{-1} (\boldsymbol{\tau}_{total} - \boldsymbol{\omega} \times (\mathbf{I}\boldsymbol{\omega})) \\
    \frac{d\mathbf{M}_{irr}}{dt} &= \text{The vector of } \frac{dM_{irr,i}}{dt} \text{ calculated in Step 3.}
\end{align}
$$
