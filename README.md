# Fully tensorial exponential PTT model

(taken from the OpenFOAM 12 User Guide - Chapter 8.3: Transport/rheology models)

$\frac{\delta \tau}{\delta t} + \nabla \bullet \left(U\tau  \right) = 2\left[\tau  \bullet \nabla U \right] - 2 \frac{\eta}{\lambda}\left(\nabla U\right) - \frac{1}{\lambda}\tau \exp \left(\frac{-\epsilon \lambda}{\eta}tr\left(\tau\right)\right)$

# One-dimensional shear simplification

## Notes

Material is sheared in x-y plane, with y being the direction of the moving plate / the telescoping motion.

## Step 1

$\frac{\delta \tau_{12}}{\delta t} + \nabla \bullet \left(U\tau_{12}  \right) = 2\left[\tau  \bullet  \dot{\gamma}\right] - 2 \frac{\eta}{\lambda}\left(\dot{\gamma} \right) - \frac{1}{\lambda}\tau$

### Special case 1: transient, no movement (quiescent melt) leading to exponential stress decay

### Special case 2: steady-state, infinitely high/long relaxation time (purely elastic edge case)

### Special case 3: steady-state, infinitely small/short relaxation time (purely viscous fluid edge case)

### Special case 4: steady-state, convective flow without shear (perfect slip)

### Special case 5: steady-state, signifcant shear, yet no convective flow across system boundaries (large diameter Couette cell) with infinitely high/long relaxation time (purely elastic edge case)

### Special case 6: the role of relaxation time on effective viscosity

# Literature
