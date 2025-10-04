# Attitude Orbit Simulator for Passive AOCS


## Equations

Angular velocity dampening by magnetic fields:

$$
\frac{d\boldsymbol{\omega}}{dt} = \mathbf{I}^{-1} \left[ \mathbf{m}_p \times \mathbf{B}(t) + \sum_{i=1}^{N} \left( \frac{V_{h,i}}{\mu_0} M_i(t) \mathbf{b}_i \right) \times \mathbf{B}(t) - \boldsymbol{\omega} \times (\mathbf{I}\boldsymbol{\omega}) \right]
$$

Simplified inertia tensor:

$$
I_x = \frac{1}{12} m (b^2 + c^2), \quad
I_y = \frac{1}{12} m (a^2 + c^2), \quad
I_z = \frac{1}{12} m (a^2 + b^2)
$$

## Dependencies

```sh
sudo geographiclib-get-magnetic wmm2020
```

```sh
sudo geographiclib-get-magnetic all
```
