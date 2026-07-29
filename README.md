## Overview
<p align="center">
  <img src="assets/demo.gif" alt="FSI simulation of vibrating robot" width="500">
</p>

This repository contains a complete OpenFOAM case for a **two-way Fluid-Structure Interaction (FSI) benchmark** based on a Vibration-driven robot (VR). The model simulates a cylinder containing an internal oscillating mass, which generates oscillations and causes the cylinder to move through a viscous incompressible fluid. This problem has a **known analytical solution**, making it an ideal test case for code verification and validation of FSI solvers.

##  Problem Description

### Physical Model

The vibrating robot consists of a rigid cylindrical body (VR) with an internal mass (IM) that oscillates along a circular path of radius `a` inside the body.

**Motion of the internal mass:**

$$
X = a\sin\Phi, \qquad Y = a\cos\Phi, \qquad \Phi(t) = \phi\cos(\omega t)
$$

where $X, Y$ are the VM displacements, $\phi$ is the oscillation amplitude, and $\omega$ is the oscillation frequency.

**Rigid-body dynamics** (6 DOF) – equations for the VR's translational velocity components $V_x, V_y$ and angular velocity $\Omega$:
**Rigid-body dynamics (3 DOF)**

$$
I_{\text{b}}\dot{\Omega} - m_{\text{im}}a^2\ddot{\Phi} = M + m_{\text{im}}(Y\dot{V}_x - X\dot{V}_y),
$$

$$
m_{\text{VR}}\dot{V}_x + m_{\text{im}}\ddot{X} = F_x,
$$

$$
m_{\text{VR}}\dot{V}_y + m_{\text{im}}\ddot{Y} = F_y,
$$

where $I_{\text{b}}$ is the moment of inertia of the body, $m_{\text{VR}}$ is the mass of the VR body, $m_{\text{im}}$ is the internal mass, and $(F_x, F_y)$ and $M$ are the hydrodynamic force components and moment from the fluid.

**Fluid dynamics** – incompressible Navier–Stokes equations for the velocity field $\mathbf{u} = (u_x, u_y)$ and pressure $p$:

$$
\frac{\partial\mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu\Delta\mathbf{u}, \qquad \nabla\cdot\mathbf{u} = 0,
$$

with $\rho$ and $\nu$ being the fluid density and kinematic viscosity.

**Hydrodynamic forces and moment** on the body are obtained by integrating the stress tensor $\boldsymbol{\sigma}$ over the surface $C$:

$$
\mathbf{F} = \oint_C \boldsymbol{\sigma}\cdot\mathbf{n} \, dl, \qquad
M = \oint_C \mathbf{r} \times (\boldsymbol{\sigma}\cdot\mathbf{n}) \, dl,
$$

where $\mathbf{n}$ is the outward normal to the surface.

This formulation provides a complete two‑way FSI coupling: the body motion influences the fluid (through moving boundaries), and the fluid exerts forces and moments that determine the body’s trajectory.

### Why This Is an Ideal Benchmark

| Feature | Status |
|:--------:|:------:|
| Bio-mimetic propulsion (flapping-wing principle) | ✓ |
| Full two-way FSI coupling | ✓ |
| Viscous incompressible flow (unsteady Navier-Stokes) | ✓ |
| Analytical solution available | ✓ |

---

##  Analytical Solution

An asymptotic solution has been derived for the steady-state motion under the assumptions:
- Small oscillation angles: $\phi \ll 1$
- Small mass ratio: $\gamma = m_{\text{im}}/m_{\text{VR}} \ll 1$

### Key Analytical Results

**Dimensionless amplitudes of translational ($\kappa$) and rotational ($\Theta$) oscillations:**

$$
\kappa = \frac{\gamma \phi}{|a_V(\beta)|}, \quad \Theta = \frac{\gamma \phi}{|a_\Omega(\beta, \alpha)|}
$$

where:

$$
a_V(\beta) \approx 2+\frac{4}{\sqrt{i\beta}}+\frac{2}{i\beta}, \quad
a_\Omega(\alpha, \beta) \approx \alpha+\frac{2}{\sqrt{i\beta}}+\frac{3}{i\beta}, \quad 
\alpha = 1 - \gamma
$$

**Steady-state (cruising) velocity:**

$$
u_{st} \approx 0.5 \frac{\phi^2 \gamma^2}{\alpha \kappa} \frac{\beta \alpha - 5.5\sqrt{\beta \alpha}}{1.33(\beta \alpha) + 12.7\sqrt{\beta \alpha} + 80.7}
$$

**Phase shifts:**

The phase shifts between the rotational ($\varphi_\Omega$) and translational ($\varphi_V$) oscillations of the robot body and the internal mass oscillation are given by:

$$
\varphi_V = \arccos\left( -\frac{\mathrm{Re}(a_V)}{|a_V|} \right), \qquad
\varphi_\Omega = \arccos\left( -\frac{\mathrm{Re}(a_\Omega)}{|a_\Omega|} \right),
$$
## Numerical Implementation in OpenFOAM
### Repository structure
```
vr-fsi-benchmark/

├── case/                        # FSI-benchmark
│   ├── 0/
│   │   ├── U
│   │   ├── p
│   │   └── params/              # initial data
│   ├── constant/
│   │   ├── polyMesh/            # baseline mesh
│   │   ├── transportProperties
│   │   └── dynamicMeshDict
│   ├── system/
│   │   ├── controlDict
│   │   ├── fvSchemes
│   │   ├── fvSolution
│   │   └── decomposeParDict
│   ├── Allrun                   # main automation script
│   └──  data.py                 # post-proccessing (in Allrun)  
├── libs/                        # Additional libraries
│   ├── externalForceSixDoF/
│   ├── cruising_auto_new/
│   └── cruising_auto/
├── meshes/                      # different meshes
│   ├── baseline/
│   ├── coarse/
│   ├── medium/
│   └── fine/
└── README.md                    
```
### Solver Configuration

| Component | Setting |
|-----------|---------|
| **Solver** | `pimpleFoam` (transient, PIMPLE algorithm) |
| **Time scheme** | Crank-Nicolson (2nd order) |
| **FSI coupling** | Two-way coupling (explicit) |
| **Motion** | Rigid body motion (6 DOF) |

### Custom class `externalForce1`

The standard OpenFOAM FSI tools do not provide a direct way to apply a local force from an internal mass. Therefore a custom class `externalForce1` was developed, inheriting from `sixDoFRigidBodyMotionRestraint`. It allows the force to be given as a function of time and automatically computes the resulting moment about the centre of mass.

**Simplified class code:**

```cpp
void externalForce1::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    scalar t = motion.time().timeOutputValue();
    restraintForce = externalForce1_().value(t);
    restraintMoment = (location_ + motion.centreOfRotation()) ^ restraintForce;
}

bool externalForce1::read(const dictionary& sDoFRBMRDict)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);
    sDoFRBMRCoeffs_.readEntry("location", location_);
    externalForce1_ = Function1<vector>::New("force", sDoFRBMRCoeffs_);
    return true;
}

class externalForce1 : public sixDoFRigidBodyMotionRestraint
{
    autoPtr<Function1<vector>> externalForce1_;
    vector location_;
    void restrain(...) const;
    bool read(const dictionary&);
};
```
Usage in constant/dynamicMeshDict – the main configuration file for the rigid‑body motion. It includes the physical parameters from 0/params and defines a cosine force profile:

```cpp
#include "./0/params"

dynamicFvMesh   dynamicMotionSolverFvMesh;
motionSolverLibs (sixDoFRigidBodyMotion sixDoFexternalForce);
motionSolver    sixDoFRigidBodyMotion;

// ... mass, inertia, constraints ...

restraints
{
    force
    {
        sixDoFRigidBodyMotionRestraint externalForce1;
        location    (1.0 0.0 0.0);
        force 
        {
            type        cosine;
            frequency   #eval{1.0/(2.0*$pi*$gamma) };
            amplitude   #eval{$mass*$phi/$gamma};
            scale       (0 1 0);
            level       (0 0 0);
        }
    }
}
```
### Integration of the equations of motion

The rigid-body dynamics equations are integrated in time using the **Newmark method** with parameters:

- `gamma* = 0.5`
- `beta* = 0.25`

This method was chosen because it provides accuracy and stability for this particular problem compared to alternative symplectic schemes.

### Mesh
The case uses a conformal mesh generated with an external tool (stored in constant/polyMesh). The mesh is refined near the cylinder to capture the boundary layer and accurately integrate forces and moments. For convergence studies, alternative meshes are provided in the meshes/ folder (coarse, medium, fine).


### Case parameters

All physical parameters are stored in `0/params`. Each file contains a single value. The following table summarises them:

| File | Symbol | Meaning | Typical value |
|------|--------|---------|---------------|
| `Vel` | `u0` | initial velocity  | set as `-u_st·κ/γ` |
| `gamma` | `γ` | mass ratio | 0.6 |
| `beta` | `β` | frequency parameter | 1000–4000 |
| `phi` | `φ` | oscillation amplitude (rad) | 0.2–0.85 |
| `pi` | `π` | constant | 3.141592654 |
| `nu1` | `ν₁` | derived `1/(βγ)` | computed automatically |

> **Important:**  
> The initial velocity must be carefully chosen. It is **not** `u_st` itself, but `-u_st · κ / γ`. This ensures consistency with the solver’s non‑dimensional form. For reference, the analytical `u_st` is given in the equation above. Vel is negative because it specifies the incoming flow velocity.
> Additionally, it should be emphasized that this initial velocity is deliberately set to differ from the steady‑state (cruise) solution. It is chosen in this way to demonstrate the robustness of the methodology for determining the cruise speed.


A typical `params` file looks like:

```cpp
Vel -0.179;    // must be -u_st * kappa / gamma
gamma 0.6;
beta 1000;
phi 0.85;
pi 3.141592654;
nu1 #calc"1.0/$beta/$gamma";
```
## Running the case

Prerequisites: OpenFOAM 2112 or newer (as it was tested on this version), Python 3 with numpy, matplotlib, pandas**

---

Clone the repository and enter the directory:
```
bash
git clone https://github.com/VDAnisimov/vr-fsi-benchmark.git
cd vr-fsi-benchmark
```
#### Additional libraries

This case requires two custom OpenFOAM libraries, which are included in the repository:

- **`externalForceSixDoF`** – implements the `externalForce1` restraint used in `dynamicMeshDict`. It applies the force from the internal oscillating mass to the rigid body at a specified point, automatically computing the resulting moment.

- **`cruising_auto`** – automatically finds the steady‑state (cruising) velocity of the vibrating robot. To do this, it uses the secant method (chord method): it iteratively adjusts the velocity until the average hydrodynamic force over one period becomes zero.
  
**Note on library versions:**
> - `cruising_auto` (original) – for OpenFOAM **2112 and older**.
> - `cruising_auto_new` – for OpenFOAM versions **newer than 2112**.

##### Compilation

Before running the case, these libraries must be compiled. Open a terminal, load the OpenFOAM environment, and execute the following commands **from the root directory of the repository** (where the `libs/` folder is located):

```bash
# Compile сruising_auto
cd libs/cruising_auto
wmake
cd ../..

# Compile externalForceSixDoF
cd libs/externalForceSixDoF
wmake
cd ../..
```
## Parallel execution

By default, the case is configured to run on **16 processors** using the `hierarchical` decomposition method with the split `(4 4 1)` in the `n` entry (settings in `case/system/decomposeParDict`). You can change the number of cores, the splitting pattern. The `Allrun` script automatically reads the `numberOfSubdomains` value and launches the solver with the appropriate number of MPI processes.

If you prefer to run in serial, simply set `numberOfSubdomains 1;` and `./Allrun` will handle the rest.
## Run

(Optional) Edit the parameters in case/0/params/ if needed. Then navigate to the case folder and run:
```
bash
./Allrun
```
### Monitoring the solution

During the simulation, you can monitor the progress of the calculation using:

- **The OpenFOAM log file**:  
The `Allrun` script writes the solver output to `case/log.pimpleFoam`.

- **The post‑processing output**:
`case/postProcessing/BCcontrol/0/coefficient.dat` contains columns:
  - first column: time,
  - middle column: aver.force ,
  - rightmost column: velocity (cruising).

These outputs help you verify that the solution is converging to the steady-state (cruising) regime.

## Numerical results summary

After the post‑processing script finishes, a file `results_table.txt` is created in the case directory. It contains the numerical values of the key characteristics extracted from the simulation:

- `beta` – frequency parameter,
- `gamma` – mass ratio,
- `phi` – oscillation amplitude,
- `u_st` – steady‑state (cruising) velocity,
- `kappa` – translational oscillation amplitude,
- `Theta` – rotational oscillation amplitude,
- `phi_V/pi` – phase shift of linear velocity (normalised by π),
- `phi_Omega/pi` – phase shift of angular velocity (normalised by π).

These numerical values can be directly compared with the analytical predictions from the **Analytical Solution** section. For the given `beta`, `gamma`, and `phi`, substitute them into the formulas for $\kappa$, $\Theta$, $u_{st}$, and the phase shifts $\varphi_V$, $\varphi_\Omega$ to obtain the reference values. The relative errors are typically within the expected ranges (≤ 3% for $\kappa$, ≤ 5% for $\Theta$, ≤ 10% for $u_{st}$).
