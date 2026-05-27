# Per-cell momentum-wind update — derivation

This document walks through the math that lives inside one cell of
the structured grid used by
`mechatree.wind._momentum_wind_kernel.compute_momentum_wind` (the
in-repo native solver) and DendroFlow's
`KEpsilonActuatorDiskModel` (the reference). It follows
Loukas Stamoulis's *"Flow through an array of cylinders at high Re"*
(`3D_model_04_05_26.pdf`, sections §1.1–§1.3) and shows that the
code is consistent with the derivation. One small typo in the PDF's
eq (19) is noted at the end.

Setup notation:

- Cubic control volume `CV` of edge `H`, frontal area `S = H²`,
  volume `V = H³`. Wind enters the upwind (`-x`) face and exits the
  downwind (`+x`) face.
- A cylinder of diameter `D` and length `L`, with axis in the `x–z`
  plane at angle `i` to `+x` (so `cos i = n_x`, `sin i = √(n_y² + n_z²)`
  for axis components `(n_x, n_y, n_z)`).
- Fluid density `ρ` is normalised to 1 in the unit-free code.
- At the six faces of `CV`, the streamwise component is `U_in`
  (upwind) and `U_out` (downwind), with lateral fluxes `V_front,
  V_back` through the `±y` faces and `W_up, W_down` through the
  `±z` faces.

---

## 1. Per-cylinder force kernel (Nat Comms 2017, vector form)

Eloy et al. (2017) write the force on a segment in **vector form**:

$$\boxed{\;\mathbf{F}_\mathrm{seg} = \tfrac{1}{2}\rho U^2 d L\, \|\mathbf{t}\times\mathbf{u}\|^2\, (\mathbf{n}\times\mathbf{t}),\qquad \mathbf{n} = \frac{\mathbf{t}\times\mathbf{u}}{\|\mathbf{t}\times\mathbf{u}\|}\;}$$

where $\mathbf{t}$ is the unit axis of the segment, $\mathbf{u}$ is
the unit wind vector, $d$ and $L$ are the diameter and length, and
the drag coefficient is taken to be 1.

### 1.1 Simplification with $\mathbf{u} = \hat{\mathbf{x}}$

With $\mathbf{t}$, $\mathbf{u}$ unit vectors, $\|\mathbf{t}\times\mathbf{u}\| = \sin i$
where $i$ is the angle between them. The direction $\mathbf{n}\times\mathbf{t}$
is the unit vector $\hat{\mathbf{m}}$ in the $(\mathbf{t},\mathbf{u})$
plane perpendicular to $\mathbf{t}$:

$$\hat{\mathbf{m}} = \frac{\mathbf{u} - (\mathbf{t}\cdot\mathbf{u})\,\mathbf{t}}{\sin i} = \frac{\mathbf{u}_\perp}{|\mathbf{u}_\perp|},\quad |\mathbf{u}_\perp| = \sin i$$

So the force becomes

$$\mathbf{F}_\mathrm{seg} = \tfrac{1}{2}\rho U^2 d L\, \sin^2 i\, \hat{\mathbf{m}} \;=\; \tfrac{1}{2}\rho U^2 d L\, \sin i\, \mathbf{u}_\perp$$

with magnitude

$$|\mathbf{F}_\mathrm{seg}| = \tfrac{1}{2}\rho U^2 d L\, \sin^2 i$$

and the **streamwise** component (projection onto $\hat{\mathbf{u}}$,
which is what feeds the per-cell momentum-wind update in §3):

$$\boxed{\;F_D = \mathbf{F}_\mathrm{seg}\cdot \hat{\mathbf{u}} = \tfrac{1}{2}\rho U^2 d L \, \sin^3 i \;}$$

(because $\hat{\mathbf{m}} \cdot \hat{\mathbf{u}} = \sin i$.)

### 1.2 What the kernel computes

For each branch the kernel
([`_momentum_wind_kernel.py:124-141`](../src/mechatree/wind/_momentum_wind_kernel.py#L124-L141)):

1. Pre-computes $\cos i = \mathbf{t}\cdot\hat{\mathbf{u}}$ (with
   $\hat{\mathbf{u}} = +\hat{\mathbf{x}}$ in the rotated wind frame,
   so $\cos i = t_x$).
2. Pre-computes $\mathbf{u}_\perp = \mathbf{u} - \cos i\,\mathbf{t}$
   (its magnitude is $\sin i$).
3. Stores three things on `MomentumWindResult`:
   - **`F_D_branch`** $= \tfrac{1}{2}\rho D L \sin^3 i \cdot U_\mathrm{loc}^2$
     — the streamwise scalar, summed per cell for the momentum-wind
     update.
   - **`F_N_branch`** $= |\mathbf{F}_\mathrm{seg}| = \tfrac{1}{2}\rho D L \sin^2 i \cdot U_\mathrm{loc}^2$
     — the normal-force magnitude (matches Taylor's $F_N$ in the
     high-Re limit).
   - **`F_vec_branch`** $= \tfrac{1}{2}\rho D L \sin i \cdot \mathbf{u}_\perp \cdot U_\mathrm{loc}^2$
     — the full 3-D Nat Comms force vector per branch, suitable for
     feeding into the tree's mechanics (option B in §7).

The drag coefficient $C_D$ is exposed as a knob (defaults to 1 per
Nat Comms) and multiplies all three.

No Reynolds-number dependence (the legacy Taylor 1952 kernel had a
$4/\sqrt{\mathrm{Re}}\,\sin^{3/2} i$ skin-friction correction that
vanishes at MechaTree's typical $\mathrm{Re} \approx 10^5$). No
lift force (the Nat Comms force is the only force we need — lift in
Taylor 1952 was a separate kernel and never affected the streamwise
wind solve anyway).

---

## 2. Conservation laws (incompressible integral form)

Loukas writes the integral Navier–Stokes equations over the `CV`:

$$\int_{\partial \mathrm{CV}} \mathbf{u}\cdot\mathbf{n}\,dS = 0 \tag{1}$$

$$\rho \int_{\partial \mathrm{CV}} \mathbf{u}(\mathbf{u}\cdot\mathbf{n})\,dS = \mathbf{F} \tag{2}$$

(advection only; viscous transport is handled separately by the
diffusion sweep — see §5 below).

Evaluating (1) for the forward-facing case `0 < i < π/2`,
`W_down, W_up > 0`:

$$U_\mathrm{in} + W_\mathrm{down} = U_\mathrm{out} + W_\mathrm{up} + |V_\mathrm{back}| + |V_\mathrm{front}| \tag{3}$$

So **mass conservation is closed by explicit lateral fluxes through
the `±y` and `±z` faces**, not "lost". Rearranging gives the key
identity used in §3:

$$U_\mathrm{in} - (V_\mathrm{front} + V_\mathrm{back}) + W_\mathrm{down} - W_\mathrm{up} \;=\; U_\mathrm{out} \tag{3'}$$

Evaluating (2) in the streamwise direction for the same case:

$$\frac{F_D}{\rho S} = U_\mathrm{in}^2 - U_\mathrm{out}^2 - U_\mathrm{in}(V_\mathrm{front} + V_\mathrm{back}) + U_\mathrm{in} W_\mathrm{down} - U_\mathrm{in} W_\mathrm{up} \tag{5}$$

The lateral momentum balances (eqs 6, 7) give the lift force
components but do not feed back into the streamwise update.

---

## 3. The momentum-wind velocity update

Factor `U_in` out of the second through fifth terms of (5):

$$\frac{F_D}{\rho S} = U_\mathrm{in}\big[U_\mathrm{in} - (V_\mathrm{front}+V_\mathrm{back}) + W_\mathrm{down} - W_\mathrm{up}\big] - U_\mathrm{out}^2$$

The bracketed expression is exactly the left-hand side of (3'), so
substitute `U_out` for it:

$$\boxed{\;\frac{F_D}{\rho S} \;=\; U_\mathrm{in} U_\mathrm{out} - U_\mathrm{out}^2 \;} \tag{18}$$

This is **what makes the per-cell update tractable**: mass
conservation lets us eliminate the lateral fluxes entirely. The
remaining equation involves only `U_in`, `U_out`, and the per-cell
drag.

(18) is a quadratic in `U_out`:

$$U_\mathrm{out}^2 - U_\mathrm{in} U_\mathrm{out} + \frac{F_D}{\rho S} = 0$$

$$U_\mathrm{out} = \tfrac{1}{2}\bigg[\,U_\mathrm{in} \pm \sqrt{U_\mathrm{in}^2 - \tfrac{4 F_D}{\rho S}}\,\bigg]
\;=\; \tfrac{U_\mathrm{in}}{2}\bigg[\,1 \pm \sqrt{1 - \tfrac{4 F_D}{\rho S \, U_\mathrm{in}^2}}\,\bigg]$$

Take the `+` root (the physical branch — `U_out > U_in/2` for
moderate drag; the `−` root is the unphysical wake-collapse branch):

$$\boxed{\;U_\mathrm{out} \;=\; \tfrac{U_\mathrm{in}}{2}\,\bigg[\,1 + \sqrt{1 - \tfrac{4 F_D}{\rho S \, U_\mathrm{in}^2}}\,\bigg]\;} \tag{19, corrected}$$

This matches the code:
[`_momentum_wind_kernel.py:171-175`](../src/mechatree/wind/_momentum_wind_kernel.py#L171-L175)

```python
denom = H**2 * U_in**2 + 1e-30      # = ρ S U_in² with ρ = 1
disc  = 1 - 4 * F_D_cell / denom
U_out = 0.5 * U_in * (1 + sqrt(disc))
```

Sanity checks:

- `F_D = 0` ⇒ discriminant = 1 ⇒ `U_out = U_in`. ✓
- `F_D → ¼ ρ S U_in²` ⇒ discriminant → 0 ⇒ `U_out → ½ U_in`. The
  momentum-wind saturation floor: the most one cell can extract is
  half the inflow.
- `F_D > ¼ ρ S U_in²` ⇒ discriminant negative. The kernel clips with
  `np.clip(disc, 0, 1)` so the cell just saturates at `½ U_in` — not
  unphysical, just the upper bound on per-cell extraction.

> **Typo flag.** The PDF's eq (19) shows the discriminant as
> $\sqrt{1 - 4 F_D / (\rho S)}$ — without the $U_\mathrm{in}^2$ in
> the denominator. That term is needed: $F_D / (\rho S)$ has units
> of velocity², so the discriminant only becomes dimensionless when
> divided by $U_\mathrm{in}^2$. The full algebra above shows it
> falls out of the quadratic solution naturally. The code has it
> right.

### 3.1 Small-drag (linear) limit

When `F_D ≪ ρ H² U_in²` — the regime where each cell extracts only a
small fraction of the momentum — expand the discriminant. Let
`ε = 4 F_D / (ρ S U_in²) ≪ 1`. Then

$$\sqrt{1 - \varepsilon} \;\approx\; 1 - \tfrac{\varepsilon}{2} - \tfrac{\varepsilon^2}{8} + \mathcal{O}(\varepsilon^3)$$

Substitute into (19, corrected):

$$U_\mathrm{out} \;\approx\; \tfrac{U_\mathrm{in}}{2}\bigg(2 - \tfrac{\varepsilon}{2}\bigg) \;=\; U_\mathrm{in} - \tfrac{U_\mathrm{in}\,\varepsilon}{4}$$

Substituting `ε` back:

$$\boxed{\;U_\mathrm{out} \;\approx\; U_\mathrm{in} \;-\; \frac{F_D}{\rho H^2 U_\mathrm{in}}\;\quad\text{(linear in }F_D\text{)}\;}$$

Equivalently, the velocity deficit `ΔU = U_in - U_out` is **linear
in the per-cell drag**:

$$\Delta U \;\approx\; \frac{F_D}{\rho H^2 U_\mathrm{in}}$$

This is the expected physical scaling: small drag → small momentum
loss → small velocity deficit, proportional to the drag and
inversely proportional to the available momentum flux
`ρ H² U_in`.

---

## 4. Equivalent obstacle for an array of cylinders (§1.3, eq 38)

If several cylinders sit in one cell, the per-cylinder drags don't
just add up — downstream cylinders sit in the wakes of upstream ones
and contribute less than they would in isolation. Loukas derives the
effective drag for an array by partitioning the cell into `N³`
sub-cells, applying eq (19) sequentially per sub-cell column,
expanding to second order in `Ψ`, and re-aggregating. The result
(his eq 38):

$$F_D^\mathrm{eq} \;=\; \sum_{\mathrm{cyl}} F_D^\mathrm{cyl}
\;-\; \frac{1}{\rho H^2 U_\infty^2}\bigg(\sum_{\mathrm{cyl}} F_D^\mathrm{cyl}\bigg)^2$$

The first term is the naive sum; the negative second term is the
**mutual-sheltering correction**.

The native kernel currently uses the naive sum
`F_D = Σ F_D^cyl` and skips the correction (it's second order in
`Ψ` per §3.1, so it vanishes in the small-drag regime and adds a
modest correction in the moderate-drag regime). Dropping it follows
the "as simple as possible" principle for the MechaTree port. If
profiling on dense-canopy cases shows it matters, add it back as one
line:

```python
F_D_eq = F_D_cell - F_D_cell**2 / (H**2 * U_in**2)
disc   = 1 - 4 * F_D_eq / (H**2 * U_in**2)
```

DendroFlow's reference kernel keeps the correction; the in-repo
port does not. Worth a benchmark before deciding.

---

## 5. Cross-stream diffusion (turbulent mixing)

After the per-cell momentum-wind update, the slice is smoothed by an
implicit 4-neighbour diffusion that represents turbulent
**cross-stream** mixing — i.e., mixing in the `(y, z)` plane at
fixed `x`, not in the streamwise direction.

$$U^\mathrm{diff}(k, j) \;=\; \frac{U(k, j) + w \cdot \big(U(k, j-1) + U(k, j+1) + U(k-1, j) + U(k+1, j)\big)}{1 + 4 w}$$

with diffusion weight

$$w \;=\; \frac{\nu_\mathrm{diff}}{|U(k, j)|\,H}$$

and edge cells clamped to themselves (`np.pad(..., mode='edge')`).

### 5.1 From k-ε to a single knob

DendroFlow's reference solver computed `ν_t` (the eddy viscosity)
from a full k-ε closure with three constants `(I_turb, l_mix, C_μ)`
plus the inflow speed `U_∞`. Substituting one step into the next:

$$\nu_t \;=\; C_\mu \cdot \frac{k_\mathrm{turb}^2}{\varepsilon}
\;,\quad
k_\mathrm{turb} = \tfrac{3}{2} U_\infty^2 I_\mathrm{turb}^2
\;,\quad
\varepsilon = \frac{C_\mu^{3/4} k_\mathrm{turb}^{3/2}}{l_\mathrm{mix}}$$

$$\Rightarrow\quad \nu_t \;=\; C_\mu^{1/4} \sqrt{\tfrac{3}{2}} \cdot U_\infty \cdot I_\mathrm{turb} \cdot l_\mathrm{mix}$$

With `C_μ = 0.09` baked in, the constant prefactor is
`α = C_μ^{1/4} √(3/2) ≈ 0.671`. So the whole k-ε apparatus
collapses to

$$\boxed{\;\nu_t \;=\; \alpha\, U_\infty\, I_\mathrm{turb}\, l_\mathrm{mix} \;\approx\; 0.671\,U_\infty I_\mathrm{turb} l_\mathrm{mix}\;}$$

— **a single number**. At the legacy defaults
`I_turb = 0.1, l_mix = 0.5, U_∞ ≈ 1`, that's `ν_t ≈ 0.034`. In
MechaTree's natural units (`L = U = 1, H = 1`), the diffusion weight
is `w = ν_t / (|U| H) ≈ 0.034`.

The native kernel exposes that one number directly as `nu_diff` —
no `I_turb`, no `l_mix`, no `C_μ`. Three magic constants replaced
with one dimensionless knob.

### 5.2 What `nu_diff` controls

`nu_diff` is the lever for the parameter sweep:

- Larger `nu_diff` ⇒ larger `w` per cell ⇒ more lateral smoothing ⇒
  **wider wake** that recovers faster downstream.
- Smaller `nu_diff` ⇒ tighter, narrower wake.
- `nu_diff = 0` disables diffusion entirely: the wake stays exactly
  the shape it had when it emerged from the canopy and gets advected
  unchanged downstream.

Try `0.01 / 0.03 / 0.1` to see how aggressively the wake spreads.

### 5.3 Is the model "just diffusion along x"?

**No.** The streamwise (`x`) direction carries advection + a drag
sink (the momentum-wind update from §3); the diffusion is purely in
the cross-stream (`y, z`) plane within each x-column. Putting (3),
(19), and (5) together as a marching PDE for `U(x, y, z)`:

$$\frac{\partial U}{\partial x} \;=\; -\underbrace{\frac{F_D^\mathrm{cell}(x, y, z)}{\rho H^2 U}}_{\text{drag sink (§3.1 linear limit)}} \;+\; \underbrace{\nabla_{yz}\!\cdot\!\big(\nu_\mathrm{diff}\,\nabla_{yz} U\big)}_{\text{lateral diffusion (§5)}}$$

— a parabolic equation with **x playing the role of "time"**
marching through the domain. There's no `∂²U/∂x²` term: no diffusion
along x.

The wake nevertheless appears to spread downstream in x because the
cross-stream diffusion is applied at every x-step and then advected
forward to the next column. So the wake's lateral width grows
monotonically with downstream distance — but the mechanism is
*lateral diffusion + x-advection*, not streamwise diffusion. Setting
`nu_diff = 0` would freeze the wake's lateral shape and carry it
downstream without spreading; the wake would still develop in x via
the drag sink, just without lateral mixing.

So: **advection (x) + sink (drag, x) + diffusion (y, z)**, not
"diffusion along x" alone.

---

## 6. Per-cell pseudocode summary

For one cell at column `i`, position `(j, k)`:

```text
# inputs: U_in[k, j], F_D_cell = Σ F_D^cyl from cylinders in this cell
denom    = ρ · H² · U_in²                   # = ρ S U_in² (ρ=1, S=H²)
disc     = clip(1 − 4 · F_D_cell / denom, 0, 1)
U_out    = ½ · U_in · (1 + √disc)           # §3 eq (19, corrected)
# then: 4-neighbour cross-stream diffusion (§5) + BCs at top z /
# lateral y → U_∞(z)
```

Per `(Nz, Ny)` slice that's one `np.bincount` (cell drag), four
lines of NumPy array arithmetic, and an `np.pad`-based diffusion —
no element-wise scalar loops, no Numba.

---

## 7. Per-branch coupling back to the tree (option B, **wired in Step 25c**)

For the canopy-aware Step-24 fixed-point loop, the wind solve and
the pruning step both need to evaluate per-branch forces:

1. **Wind solve** (this kernel) computes per-branch `U_loc` (the
   inflow at the branch's cell) and per-branch `F_N` (the normal
   force from the high-Re Nat Commun kernel above), then aggregates
   `F_D = F_N sin i` per cell for the momentum-wind update.
2. **Pruning** ([pruning.cpp](../src/mechatree/_core/pruning.cpp))
   recomputes per-branch forces from a single canopy-mean wind and
   feeds them to the Weibull failure law.

The recomputation in (2) is inconsistent with (1) — the wind solve
already knows each branch's local `U_loc` and produced its
`F_vec`/`F_N`, so collapsing to a canopy-mean throws away the spatial
structure (a sheltered understory twig and an exposed crown tip get
the same wind). **Option B** plumbs the per-branch force through:

- The kernel exposes `MomentumWindResult.F_vec_branch` (the full
  per-branch force vector) and `U_branch` (per-branch local wind
  magnitude).
- `Branch` carries two dedicated fields `segment_force_` and
  `segment_wind_` ([branch.h](../src/mechatree/_core/branch.h)),
  written via `PyTree.set_segment_forces_batch` /
  `set_segment_winds_batch`. These are **separate** from the `force_`
  slot, which `prune` reuses as its leaves-to-trunk aggregation
  accumulator — so the pre-stored CFD force survives the walk.
- A new C++ entry `prune_with_stored_forces(tree, leaf_drag_S0,
  cauchy)` ([pruning.cpp](../src/mechatree/_core/pruning.cpp)) reads
  `segment_force_` as each branch's own woody drag (skipping the
  `wind_force` recompute) and `segment_wind_` for the per-branch
  leaf-cluster drag term.
- [`MomentumWindBridge`](../src/mechatree/wind/momentum.py)
  writes both arrays (rotated to the world frame) onto each tree after
  the solve; the prune loops dispatch to `prune_with_stored_forces`.

As of Step 26a this is the **only** momentum behaviour — there is no
canopy-mean wind mode (a 3-D solve averaged to one scalar is just a
constant wind). The canopy-mean survives only as the ε convergence
thermometer for the Step-24 fixed-point loop. The win is physics
fidelity, not speed: per-branch forces prune more than a canopy-mean
because exposed branches feel the full storm.

---

## Summary of consistency check vs. Loukas's PDF

| § PDF | What it says | What the code does | Status |
|-------|--------------|--------------------|--------|
| §1.1, eq 12–15 | Per-cylinder force kernels (Taylor 1952) | High-Re limit: `F_D = ½ C_D D L sin³ i · U²`; no lift, no Re skin-friction | ✓ Matches §1.1 limit |
| §1.1, eqs 1–7 | Mass + x-momentum integral balance over `CV` | Used to derive eq 18 below | ✓ Matches |
| §1.2, eq 18 | `F_D/(ρS) = U_in U_out - U_out²` | Implicit in the quadratic | ✓ Matches |
| §1.2, eq 19 | `U_out = ½ U_in [1+√(1 - 4F_D/(ρS))]` | `U_out = ½ U_in [1+√(1 - 4 F_D / (ρ S U_in²))]` | ⚠️ **PDF typo**: missing `U_in²` in denominator. Code uses the correct (dimensionally consistent) form |
| §1.3, eq 38 | `F_D^eq = Σ F_D^cyl - (Σ F_D^cyl)² / (ρ H² U_∞²)` | Native port **drops** the correction (second-order); DendroFlow reference keeps it | ⚠️ Intentional simplification — see §4 |
| (not in PDF) | Turbulent diffusion | k-ε closure collapses to single `nu_diff` — see §5 | OK; documented |
