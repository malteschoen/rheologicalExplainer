# Fully tensorial exponential PTT model

(taken from the OpenFOAM 12 User Guide - Chapter 8.3: Transport/rheology models)

$\frac{\delta \tau}{\delta t} + \nabla \bullet \left(U\tau  \right) = 2\left[\tau  \bullet \nabla U \right] - 2 \frac{\eta}{\lambda}\left(\nabla U\right) - \frac{1}{\lambda}\tau \exp \left(\frac{-\epsilon \lambda}{\eta}tr\left(\tau\right)\right)$

# One-dimensional shear simplification

### Step 1

$\frac{\delta \tau_{12}}{\delta t} + \nabla \bullet \left(U\tau_{12}  \right) = 2\left[\tau  \bullet  \dot{\gamma}\right] - 2 \frac{\eta}{\lambda}\left(\dot{\gamma} \right) - \frac{1}{\lambda}\tau$

# Literature
