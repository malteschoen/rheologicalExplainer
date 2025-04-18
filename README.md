# Fully tensorial exponential PTT model

(taken from the OpenFOAM 12 User Guide - Chapter 8.3: Transport/rheology models)

$\frac{\delta \tau}{\delta t} + \nabla \bullet \left(U\tau  \right) = 2\left[\tau  \bullet \nabla U \right] - 2 \frac{\eta_0}{\lambda}\left(\nabla U\right) - \frac{1}{\lambda}\tau \exp \left(\frac{-\epsilon \lambda}{\eta_0}tr\left(\tau\right)\right)$

## Notes:


# One-dimensional (scalar) shear simplification neglecting first normal stress difference

## Notes

Material is sheared in x-y plane, with y being the direction of the moving plate / the telescoping motion.
We choose to ignore all stresses other than $\tau_{xy}$, hence  $tr\left(\tau\right)=0$, leading to the expoential term becoming  $1$

## First transformation: from $\nabla U$ to $\dot{\gamma}_{xy}$

```math
\frac{\delta \tau_{xy}}{\delta t} + \nabla \bullet \left(u_y  \tau_{xy}  \right) = 2\left[\tau_{xy}  \bullet  \dot{\gamma}_{xy}\right] - 2 \frac{\eta_0}{\lambda}\left(\dot{\gamma}_{xy} \right) - \frac{1}{\lambda}\tau_{xy}
```

## Second transformation: usage of inner and outer products
## Third transformation: from $\nabla$ to derivatives

### Special case 1: transient, no movement (quiescent melt) leading to exponential stress decay

### Special case 2: steady-state, infinitely high/long relaxation time (purely elastic edge case)

### Special case 3: steady-state, infinitely small/short relaxation time (purely viscous fluid edge case)

### Special case 4: steady-state, convective flow without shear (perfect slip)

### Special case 5: steady-state, signifcant shear, yet no convective flow across system boundaries (large diameter Couette cell) with infinitely high/long relaxation time (purely elastic edge case)

### Special case 6: the role of relaxation time on effective viscosity

# Literature
