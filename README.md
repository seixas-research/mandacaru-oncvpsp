# Mandacaru Optimized Norm-Conserved Vanderbilt Psedopotentials (ONCVPSP)

Optimized norm-conserving Vanderbilt pseudopotentials (ONCVPSP) for
[Mandacaru](https://github.com/seixas-research/mandacaru): one file per
element for every element with **Z ≤ 92** (H through U). Mandacaru generates
them from scratch with its own all-electron radial solver and its
`mandacaru.pseudopotentials.oncv` module; nothing here comes from another code.

The datasets live here, not in the Mandacaru package, because of their size
(about 125 MB).

## Contents

| Folder | Datasets | Reference atom |
|---|---|---|
| `lda/` | 92 (H–U) | LDA, scalar-relativistic, nonlinear core correction |

Mandacaru reads `<checkout>/<xc>/`, so a PBE set would go in `pbe/` beside
`lda/`.

## Using the datasets

```bash
git clone https://github.com/seixas-research/mandacaru-oncvpsp
mandacaru --set-oncvpsp /path/to/mandacaru-oncvpsp   # writes MANDACARU_ONCVPSP_PATH
# open a new terminal, then
mandacaru --pseudo-status
```

The family is selected as a basis, `basis="ONCVPSP"` (alias `"ONCV"`) or
`basis={"name": "ONCVPSP", "size": "DZP"}`:

```python
from ase.build import molecule
from mandacaru import Mandacaru

atoms = molecule("H2O")
atoms.center(vacuum=4.0)
atoms.calc = Mandacaru(method="adapt-vqe", basis="ONCVPSP", h=0.25)
atoms.get_potential_energy()
```

Without `MANDACARU_ONCVPSP_PATH` an ONCVPSP calculation stops before it
starts, with a `LibraryPathError` that names the command above. The basis
option `directory="..."` points a single run at another folder of datasets.

## Construction

Following D. R. Hamann, *Phys. Rev. B* **88**, 085117 (2013):

1. a scalar-relativistic LDA atom; two partial waves per `l` (the bound level
   and one scattering energy above it);
2. smooth pseudo waves as spherical-Bessel expansions inside `r_c`, with
   generalized norm conservation and the residual kinetic energy beyond
   `q_c = 5 Bohr⁻¹` minimized exactly;
3. a polynomial local potential inside `r_cl`, which follows the **largest**
   cutoff; each channel's two projectors reach out to `r_cl`, where they are
   `(V_AE − V_loc) φ`, so no channel sits in the bare all-electron well;
4. a 2×2 coupling matrix per `l`, and unscreening with a partial core density
   (nonlinear core correction, Louie, Froyen and Cohen 1982).

## Checks, and the flagged elements

Every dataset was checked when it was generated:

- **ghost states:** the spectrum of every channel, with and without
  projectors, against the all-electron reference;
- **scattering:** the phase `arctan L(E)` of the logarithmic derivative,
  compared with the all-electron atom at the projector radius, within
  0.05 rad over ε ± 0.5 Ha and 0.3 rad over ε ± 1 Ha.

When the default construction failed a check, the generator tried raised
local potentials and balanced cutoffs. **No element holds a ghost state.** An
element that no repair cleans is still written, with its defect recorded in
the file; loading it raises a `GhostStateWarning`, and its `repr` says
`SCATTERING OFF`.

| Elements | Defect |
|---|---|
| Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy | `f`-channel phase off by 0.06–0.11 rad |
| Tl, Pb, Bi, Po, At, Rn | phase of the semicore 4f channel off by 0.49–0.65 rad |

For the compact 4f channel this phase is measured at a radius where the 4f
wave has all but vanished, so the second row may overstate the error.

## Generating the datasets

With a Mandacaru development install and `MANDACARU_ONCVPSP_PATH` set:

```bash
mandacaru-build --pp ONCV --relativistic --xc LDA --all --workers 5 --check --ghosts flag --install
```

`--install` writes into `$MANDACARU_ONCVPSP_PATH/lda/`. The radial kernels
run in C; the full set takes under two hours on 5 cores.

## File format

Each `<Symbol>.parquet` is a self-describing Mandacaru pseudopotential record
(format `mandacaru-pseudopotential`, version 2, family `oncvpsp`) on a
3000-point radial grid (0.01 Bohr out to 30 Bohr). It holds the pseudo waves,
the projectors and coupling matrices, the local potential (screened and
unscreened), the valence and partial core densities, plus the generation
record (functional, relativity, core correction) and any recorded defects.
All quantities are in atomic units. Read one with
`mandacaru.pseudopotentials.io.load_pseudopotential(path)`.

## License

MIT, see [LICENSE](LICENSE).
