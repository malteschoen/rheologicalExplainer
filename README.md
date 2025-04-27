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

# Stabilizations according to rheoTool user guide and others

| Type | LHS: Divergence of polymeric stress contribution | RHS: Divergence of polymeric stress contribution - overbrace indicates linear interpolation from cell centre instead of just 'dropping in' cell centre value | RHS: Widehat indicates 'special second-order derivative'. $\eta_p$ is a habitual scaling factor, but you might as well put in anything else.   | RHS: Divergence of solvent contribution (from solvent viscosity and deformation gradient) plus stabilization terms |
| ---- | ---- | ---- | ---- | ---- |
| none (rheoTool)     | $\nabla \cdot \tau =$ | $\nabla \cdot \overbrace{\tau}$  | ---- |$\nabla \cdot (\eta_s \nabla U)$ |
| BSD (rheoTool)      | $\nabla \cdot \tau =$ | $\nabla \cdot \overbrace{\tau}$  | $-\nabla \cdot (\eta_p \nabla U)$  | $\nabla \cdot  (\eta_s + \eta_p) \nabla U$ |
| coupling (rheoTool) | $\nabla \cdot \tau =$ | $\nabla \cdot \overbrace{\tau}$  | $-\widehat{\nabla \cdot (\eta_p \nabla U)}$ | $\nabla \cdot  (\eta_s + \eta_p) \nabla U$ |

- Do note that for a suffciently fine mesh $\tau$ and $\overbrace{\tau}$ are equal!
- Do note that for a suffciently fine mesh the other added terms also cancel out.
- Do note that these formulae give NO information about what the value of $\tau$ would be. That is being handled by the constitutive equation. The equations above rather add a bit of 'skillful noise' to help with stability.


# Stabilization (?) as implemented in the new momentumTransport models of OpenFOAM
### Unadulterated source code 
```
template<class BasicTurbulenceModel>
tmp<Foam::fvVectorMatrix>
Maxwell<BasicTurbulenceModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    return
    (!play https://www.youtube.com/watch?v=8Q-b3bLQ3jc
        fvc::div
        (
            this->alpha_*rho*this->nuM_*fvc::grad(U)
        )
      + fvc::div(this->alpha_*rho*sigma_)
      - fvc::div(this->alpha_*rho*this->nu()*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*rho*nu0(), U)
    );
}
```
### Simplifed source code
 
```
 fvc::div
        (
            nuP*fvc::grad(U)
        )
      + fvc::div(tau)
      - fvc::div(nuS*dev2(T(fvc::grad(U))))
      - fvm::laplacian(nuP+nuS, U)
```
- we assumed $\alpha = \rho= 1$
- we renamed $nuM = nuP$
- we renamed $sigma= tau$
- we exploited $nu0 = nuS+nuP$ (see definition in Maxwell.H)
- remember: fvc is explicit, fvm is implicit

| Type | LHS: Divergence of polymeric stress contribution | RHS: Divergence of polymeric stress contribution - overbrace indicates linear interpolation from cell centre instead of just 'dropping in' cell centre value | RHS: Widehat indicates 'special second-order derivative'. $\eta_p$ is a habitual scaling factor, but you might as well put in anything else.   | RHS: Divergence of solvent contribution (from solvent viscosity and deformation gradient) plus stabilization terms |
| ---------- | ---- | ---- | ---- | ---- |
| as code    | fvc::div(tau)| - fvc::div(nuS*dev2(T(fvc::grad(U)))) |fvc::div(nuP*fvc::grad(U)) | -fvm::laplacian(nuP+nuS, U)|
| BSD from above as rference |  $\nabla \cdot \tau =$ | $\nabla \cdot \overbrace{\tau}$  | $-\nabla \cdot (\eta_p \nabla U)$  | $\nabla \cdot  (\eta_s + \eta_p) \nabla U$ |

# Literature
