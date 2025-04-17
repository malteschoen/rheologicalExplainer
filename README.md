# rheologicalExplainer

this is an exponential PTT model (taken from the OpenFOAM 12 User Guide - Chapter 8.3: Transport/rheology models)

$ \frac{\delta \tau}{\delta t} + \nabla \bullet \left(U\tau  \right) = 2\left[\tau  \bullet \nabla U \right] - 2 \frac{\eta}{\lambda}\left(\nabla U\right) - \frac{1}{\lambda}\tau \exp \left(\frac{-\epsilon \lambda}{\eta}tr\left(\tau\right)\right)$
