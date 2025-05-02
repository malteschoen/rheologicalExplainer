# Chapter 1: Fully tensorial exponential PTT model

(taken from the OpenFOAM 12 User Guide - [Chapter 8.3: transport and rheology models](https://doc.cfd.direct/openfoam/user-guide-v12/transport-rheology))

$\frac{\delta \tau}{\delta t} + \nabla \bullet \left(U\tau  \right) = 2\left[\tau  \bullet \nabla U \right] - 2 \frac{\eta_0}{\lambda}\left(\nabla U\right) - \frac{1}{\lambda}\tau \exp \left(\frac{-\epsilon \lambda}{\eta_0}tr\left(\tau\right)\right)$

### Notes on variables
- $\tau$ is a (symmetric) 3x3 matrix/tensor, holding normal stresses on the diagonal and shear stresses below/above it. SI unit: Pascal
- $U$ is a vector. SI unit: meter per second
- $\lambda$ is the scalar (characteristic) relaxation time. SI unit: seconds
- $\eta_0$ is the scalar zero-shear viscosity. SI unit: Pascal-seconds
- $\epsilon$ is the scalar "elongational parameter". Dimensionless.

### Notes on mathematical operators
- "$\nabla U$" is the gradient of the velocity vector; i.e. a 3x3 tensor/matrix usually known as gradU
- More generally speaking, $\nabla$ increases the rank r of its victim by one (scalars become vectors, vectors become matrices). It is used to look at the local gradient.
-
- "$\bullet $" indicates the inner product (takes two operands of identical rank r and returns an output that also has the rank r). 
- Familiar examples would be simple scalar multiplication or XYZ
-
- "$\nabla \bullet$" indicates the divergence operator which decreases the rank of its victim by one (matrices become vectors, vectors become scalars). 
- It is used to look at the rate of destruction/generation - a divergence-less velocity field (i.e. incompressible flow) has neither generation nor destruction of flowrate inside it. 
-
- $tr$ is the trace operator. Only applicable to matrices, it returns the scalar sum of the diagonal entries.





# One-dimensional (scalar) shear simplification neglecting first normal stress difference

### Notes

Material is sheared in x-y plane, with y being the direction of the moving plate / the telescoping motion.
We choose to ignore all stresses other than $\tau_{xy}$, hence  $tr\left(\tau\right)=0$, leading to the expoential term becoming  $1$

### First transformation: from $\nabla U$ to $\dot{\gamma}_{xy}$

$\frac{\delta \tau_{xy}}{\delta t} + \nabla \bullet \left(u_y  \tau_{xy}  \right) = 2\left[\tau_{xy}  \bullet  \dot{\gamma}_{xy}\right] - 2 \frac{\eta_0}{\lambda}\left(\dot{\gamma}_{xy} \right) - \frac{1}{\lambda}\tau_{xy}$

### Second transformation: usage of divergence and inner products

$\frac{\delta \tau_{xy}}{\delta t} + \frac {\delta (u_y  \tau_{xy})}{\delta y } = 2\left[\tau_{xy}  \dot{\gamma}_{xy}\right] - 2 \frac{\eta_0}{\lambda}\left(\dot{\gamma}_{xy} \right) - \frac{1}{\lambda}\tau_{xy}$

## List of special cases: 

### Special case 1: transient, no movement (quiescent melt) leading to exponential stress decay
By eliminating all terms with $\dot{\gamma}_{xy}$ and $u_{y}$ in it, we arrive at:

$\frac{\delta \tau_{xy}}{\delta t}= \frac{-1}{\lambda}\tau_{xy}$

The solution of this differential equation must have the shape of $\tau_{xy}(t)= \exp( \frac{-t}{\lambda})$. We can conclude that the right-hand-side term $\frac{-1}{\lambda}\tau_{xy}$ represents the stress relaxation.

### Special case 2a: steady-state, infinitely high/long relaxation time (purely elastic edge case)
By setting $\lambda$ to infinity and neglecting $\frac{\delta \tau_{xy}}{\delta t}$

### Special case 3: steady-state, infinitely small/short relaxation time (purely viscous fluid edge case)
By  multiplying everything with $\lambda$ and then setting $\lambda$ to 0 , we arrive at at $\eta_0 \dot{\gamma}_{xy}=\tau_{xy}$.
This is the definition of a Newtonian fluid. 


### Special case 4: steady-state, convective flow without shear (perfect slip)
If flow through our steady-state system would occur without any shear, we would get:

### Special case 5: steady-state, signifcant shear, yet no convective flow across system boundaries (large diameter Couette cell)

### Special case 6: the role of relaxation time on effective viscosity

# Chapter 2a: Stabilizations according to rheoTool user guide and others

| Type | LHS:| RHS 1: overbrace indicates linear interpolation from cell centre instead of just 'dropping in' cell centre value | RHS 2: Widehat indicates 'special second-order derivative'. $\eta_p$ is a habitual scaling factor, but you might as well put in anything else.   | RHS 3: Divergence of solvent contribution (from solvent viscosity and deformation gradient) plus stabilization terms (?) | Helpful link |
| ---- | ---- | ---- | ---- | ---- | ---- |
| none (rheoTool) in OpenFOAM | ---- | fvc::div(tau) | ---- |fvm::laplacian(etaS,U)| ---- |
| none (rheoTool) in equations     | $\nabla \cdot \tau =$ | $\nabla \cdot \overbrace{\tau}$  | ---- |$\nabla \cdot (\eta_s \nabla U)$ | ---- |
| ---- | ---- | ---- | ---- | ---- | ---- |
| BSD (rheoTool) in OpenFOAM | ---- | fvc::div(tau)  |  -fvc::laplacian(etaP, U) | fvm::laplacian((etaP + etaS),U) | ---- |
| BSD (rheoTool) in equations     | $\nabla \cdot \tau =$ | $\nabla \cdot \overbrace{\tau}$  | $-\nabla \cdot (\eta_p \nabla U)$  | $\nabla \cdot  (\eta_s + \eta_p) \nabla U$ | ---- |
| ---- | ---- | ---- | ---- | ---- | ---- |
| coupling (rheoTool) in OpenFOAM | ---- | fvc::div(tau) | -etaP*fvc::div(fvc::grad(U)) | fvm::laplacian((etaP + etaS),U)| ---- |
| coupling (rheoTool) in equations | $\nabla \cdot \tau =$ | $\nabla \cdot \overbrace{\tau}$  | $-\widehat{\nabla \cdot (\eta_p \nabla U)}$ | $\nabla \cdot  (\eta_s + \eta_p) \nabla U$ | ---- |
| ---- | ---- | ---- | ---- | ---- | ---- |
| DEVSS (Favero) in OpenFOAM | ---| fvc::div(tau)|  -fvc::laplacian(etaP,U)| fvm::laplacian((etaP + etaS),U) |[here](https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-4.0/blob/268bb07d15d8d2de5df531f7702df54da05f15ad/src/transportModels/viscoelastic/viscoelasticLaws/Giesekus/Giesekus.C)|
| DEVSS (Favero) in equations | $\nabla \cdot \tau =$ | $\nabla \cdot \overbrace{\tau}$  | $-\nabla \cdot (\eta_p \nabla U)$ | $\nabla \cdot  (\eta_s + \eta_p) \nabla U$ |[here](https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-4.0/blob/268bb07d15d8d2de5df531f7702df54da05f15ad/src/transportModels/viscoelastic/viscoelasticLaws/Giesekus/Giesekus.C)|

- Do note that for a suffciently fine mesh $\tau$ and $\overbrace{\tau}$ are equal! What does that even mean?
- Do note that for a suffciently fine mesh the other added terms also cancel out. Which terms?
- Do note that these formulae give NO information about what the value of $\tau$ would be. That is being handled by the constitutive equation. The equations above rather add a bit of 'skillful noise' to help with stability.

### Stabilization (?) as implemented in the new momentumTransport models of OpenFOAM
#### Unadulterated source code from OF2312 [here](https://www.openfoam.com/documentation/guides/latest/api/Maxwell_8C.html)
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
    (
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
#### Simplifed source code from OF v2312
 
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
- we renamed $sigma = tau$
- we exploited $nu0 = nuS+nuP$ (see definition in Maxwell.H)
- remember: fvc is explicit, fvm is implicit

#### for reference: simplified version of laminarModel.C from OpenFOAM v2312 (see [here](https://www.openfoam.com/documentation/guides/latest/api/linearViscousStress_8C_source.html))

```
      - fvc::div(nuEff*dev2(T(fvc::grad(U))))
      - fvm::laplacian(nuEff, U)
 ```
# Chapter 2b:  Understanding OpenFOAM additions

| Type | LHS: | RHS 1: | RHS 2:  | RHS 3:  | Helpful link |
| ---------- | ---- | ---- | ---- | ---- | ---- |
| GNF (rheoTool) | ---- | fvc::grad(U) & fvc::grad(eta) |----|fvm::laplacian(eta,U)|----|
| VE, no stabil. (rheoTool) | fvc::div(tau) | ---- | ---- |fvm::laplacian(etaS,U)| ---- |
| VE, BSD stabil. (rheoTool)  | fvc::div(tau) | ---- |  -fvc::laplacian(etaP, U) | fvm::laplacian((etaP + etaS),U) | ---- |
| VE, coupling stabil. (rheoTool)  | fvc::div(tau) | ---- |  - etaP*fvc::div(fvc::grad(U)) | fvm::laplacian((etaP + etaS),U) | ---- |
| laminarModel.C (OpenFOAM v2312)      |----- | - fvc::div(nuEff* dev2(T(fvc::grad(U)))) | ----| -fvm::laplacian(nuEff, U) | ---- |
| Maxwell.C  (OpenFOAM v2312)          | fvc::div(tau) | - fvc::div(nuS* dev2(T(fvc::grad(U)))) | fvc::div(nuP*fvc::grad(U)) | -fvm::laplacian(nuP+nuS, U) | ---- |

# Chapter 3: The log-conformation approach

# Literature
