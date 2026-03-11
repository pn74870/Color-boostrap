# Color-boostrap

Numerical S-matrix bootstrap tooling for **adjoint scalar scattering with color structure**, extending the twice-subtracted dispersion framework of [Caron-Huot, Van Duong, Vanichchapongjaroen, Vichi (2021)](https://arxiv.org/pdf/2102.08951) to non-abelian representation channels.

This repository implements the color-improved dispersion-relation workflow used to extract bounds on EFT couplings from:
- crossing,
- unitarity/positivity,
- finite/large impact-parameter constraints,
- and (optionally) the graviton pole contribution in the singlet channel.

## What is in this repository

- `disp_rels_color_general.nb`  
  Main Mathematica notebook deriving/assembling the color dispersion relations and running bound computations.
- `full_bounds_2sub_int_fwd.wl`  
  Core Wolfram Language package/script with the bound pipeline (`calculateBounds`) and helper routines.
- `SDPB.m`  
  SDPB interface utilities used by the bootstrap workflow.
- `generalColorJ30m10.wdx`  
  Cached data artifact (example precomputed data).

## Physics setup (high-level)

The code uses a twice-subtracted fixed-`t` dispersion relation for
\[
\phi^a\phi^b \to \phi^c\phi^d,
\]
with explicit color projectors and crossing maps between tensor bases. Compared to the baseline 2102.08951 setup, this project:

1. keeps full representation-space structure (via matrices such as `M^+`, `W`, and `X=M^+W`),
2. isolates the improved combination where only the `t`-channel singlet graviton pole survives,
3. eliminates higher EFT Taylor data (`n\ge 3` in `s`-derivatives) in favor of dispersive moments,
4. implements the resulting improved kernels (including even/odd rep channels) numerically.

The notebook and package embody the component/vector relations and improved singlet equation described in your derivation notes.

## Prerequisites

- **Wolfram Mathematica** (or Wolfram Kernel compatible with `.wl` execution).
- **Docker** (required by the default SDPB runner).
- Access to image `wlandry/sdpb:2.5.1`.
- Sufficient CPU/RAM for high precision semidefinite programs.

The script imports `SDPB.m` and calls SDPB through Docker + `mpirun`.

## Quick start

1. Open Mathematica in this repository directory.
2. Evaluate (or load) `full_bounds_2sub_int_fwd.wl`.
3. Define your theory-specific objects from the notebook, including:
   - `dispRelHighE`, `dispLowE`,
   - improved kernels: `CimpEvenR`, `CimpOddR`, `CimpLowE`,
   - color maps: `Mp`, `W` (or model-specific variants like `MplusSO2`, `WSO2`),
   - EFT scan rules `eftVals`.
4. Call `calculateBounds[...]` (details below).

---

## How to call `calculateBounds`

`calculateBounds` is defined in `full_bounds_2sub_int_fwd.wl` and documented in comments directly above the function.

### Signature

```wl
calculateBounds[
  eftVals_,
  dispRelHighE_, dispLowE_,
  CimpEvenR_, CimpOddR_, CimpLowE_,
  fileFcVals_, integralsFileName_, fwdFileName_,
  nMax_Integer, Jmax_Integer, nFwdMax_Integer,
  Mp_, W_,
  kMax_Integer, d_Integer,
  xs_, \[Delta]b_, Bmax_, mmax_,
  polsFileName_, savePolsFile_, boundsFileName_,
  hasGrav_:True,
  excludeReps_:{},
  inclLargeImpPar_:True
]
```

### Argument guide (practical)

1. `eftVals`  
   List of replacement-rule sets. Each set defines fixed EFT inputs and leaves `g` as objective.
2. `dispRelHighE`, `dispLowE`  
   Non-improved high/low-energy dispersion relations.
3. `CimpEvenR`, `CimpOddR`  
   Improved high-energy kernels for symmetric/antisymmetric reps.
4. `CimpLowE`  
   Improved low-energy relation in terms of `B[s,t]`, `WB[s,t]=W.B[s,t]`, and graviton pole term.
5. `fileFcVals`  
   Cache file for integral helper values.
6. `integralsFileName`  
   Cache/output for computed high-energy and finite-impact-parameter integrals.
7. `fwdFileName`  
   Cache/output for forward-limit equations.
8. `nMax`, `Jmax`  
   Integral truncations (smearing-power and maximum spin).
9. `nFwdMax`  
   Maximum number of forward-limit equations included.
10. `Mp`, `W`  
   Crossing/color matrices (`abcd->adbc` and `abcd->abdc`, with `W` typically diagonal ±1).
11. `kMax`  
   Maximum Taylor order for forward-limit equation generation.
12. `d`  
   Spacetime dimension.
13. `xs`  
   Sampling points in `x` with `\[Mu]=1/(1-x)`.
14. `\[Delta]b`, `Bmax`  
   Step size and cutoff for finite impact-parameter constraints.
15. `mmax`  
   Number of subleading terms for large impact-parameter asymptotics.
16. `polsFileName`, `savePolsFile`  
   Cache/output and reuse flag for SDPB polynomial data.
17. `boundsFileName`  
   Output file for min/max bounds and optimal functionals.
18. `hasGrav` (optional, default `True`)  
   Include singlet graviton pole logic; when `False`, forward equations include the t-singlet treatment accordingly.
19. `excludeReps` (optional)  
   List of representation indices excluded from positivity matrices.
20. `inclLargeImpPar` (optional, default `True`)  
   Toggle large impact-parameter constraints.


### `getCimpLowE` / `CimpLowE` (t-channel projector-basis check)

The notebook’s improved low-energy object is the `CimpLowE` expression (this is the object used by the helper often referred to as `getCimpLowE` in analysis discussions):

```wl
CimpLowE =
  (8 G Pi Pt[1])/t
  + (1/2) (B^(0,2)[0,0] - WB^(0,2)[0,0])
  - B^(1,1)[0,0]
  - (1/2) t B^(1,2)[0,0];
```

Interpretation in the fixed `t`-channel projector basis:
- `B` is the vector of low-energy coefficients in representation space,
- `WB` is `W.B`, i.e. the crossed/sign-flipped rep vector,
- `Pt[1]` selects the `t`-channel singlet projector, so `(8 G Pi Pt[1])/t` is the explicit graviton-pole piece,
- the derivative terms are the finite EFT subtraction data retained after improved elimination of higher-`s` tails.

In other words, `CimpLowE` is the explicit low-energy side of the improved dispersion relation after projecting to the fixed `P^{adbc}` basis and truncating the subtraction sector to the derivatives needed by the improved construction.

To inspect the **explicit** form in the fixed `t`-projector basis, use:

```wl
getCimpLowE[Mp_,W_,CimpLowE_]:=Module[{symmReps,basis,wbasis,ss,tt},
  symmReps=Flatten@Position[Diagonal[W],_?Positive];
  basis=D[getB[Mp,3,symmReps],{Array[Pt,Length[Mp]],1}];
  wbasis=W . basis;
  ss=Unique["s"];
  tt=Unique["t"];
  CimpLowE/. {
    B[x_,y_]:>(basis/. {s->x,t->y}),
    WB[x_,y_]:>(wbasis/. {s->x,t->y}),
    Derivative[m_,n_][B][x_,y_]:>(D[basis/. {s->ss,t->tt},{ss,m},{tt,n}]/. {ss->x,tt->y}),
    Derivative[m_,n_][WB][x_,y_]:>(D[wbasis/. {s->ss,t->tt},{ss,m},{tt,n}]/. {ss->x,tt->y}),
    Pt[1]->D[Pt[1],{Array[Pt,Length[Mp]],1}]
  }
]
```

What this helper does:
- builds the low-energy basis vector (`basis`) from `getB[Mp,3,symmReps]`,
- builds crossed basis coefficients (`wbasis = W . basis`),
- replaces symbolic `B`, `WB`, and all needed derivatives by explicit basis-vector expressions,
- replaces `Pt[1]` by its coefficient-vector form,
- returns the concrete rep-space vector expression for `CimpLowE` in the `t`-projector basis.

Pipeline usage details:
- `calculateBounds[...]` passes `CimpLowE` into `GetLowEintegrals[CimpLowE, Mp, W, 3, nMax]` (only `k=3` is needed in this improved relation).
- `GetLowEintegrals` performs the same type of substitutions internally to assemble low-energy integral constraints before SDPB.

This is exactly the implementation-level bridge between the symbolic improved equation and the linear constraints fed to SDPB.

### Example call pattern

From notebook comments/cells, a typical invocation is:

```wl
calculateBounds[
  {
    {
      G -> (1/8)/Pi,
      g[2,1] -> 1/2,
      g[2,2] -> -(g+1)/2
    },
    {
      G -> (1/8)/Pi,
      g[2,1] -> 20,
      g[2,2] -> -(g+1)/2
    }
  },
  dispRelHighE, dispLowE,
  CimpEvenR, CimpOddR, CimpLowE,
  fileFcVals, integralsFileName, "fwdSo2.wdx",
  10, 30, 0,
  MplusSO2, WSO2,
  10, 5,
  xs, \[Delta]b, Bmax, mmax,
  "pols_so2.wdx", False,
  "bounds.wdx",
  True, {}, True
];
```

Another cell shows an SO(5)-style call with `MpSonAdj /. n->5` and `WSonAdj`.

## Caching and outputs

The workflow is cache-first:
- If `integralsFileName` exists, integrals are imported instead of recomputed.
- If `fwdFileName` exists, forward-limit equations are imported.
- If `polsFileName` exists and `savePolsFile=True`, SDPB polynomial matrices are imported.

Main result files are `.wdx` exports containing bounds and associated functional data.

## Notes on execution

- SDPB runs are launched through Docker containers from Wolfram code.
- The script uses high precision and can be computationally heavy.
- Tune `Jmax`, `nMax`, `nFwdMax`, sampling `xs`, and impact-parameter settings to trade cost vs. tightness of bounds.

## Citation

If this code contributes to your work, please cite:
- S. Caron-Huot et al., *Sharp Boundaries for the Swampland* (2021), [arXiv:2102.08951](https://arxiv.org/pdf/2102.08951),
- and describe that your implementation uses a color-extended improved dispersion framework for adjoint external scalars.
