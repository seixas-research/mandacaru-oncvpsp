# mandacaru-oncvpsp

Optimized norm-conserving Vanderbilt pseudopotentials (ONCVPSP) for
[Mandacaru](https://github.com/seixas-research/mandacaru), one file per element
for every element with **Z ≤ 92** (H through U). The datasets are generated
from scratch by Mandacaru's own LDA radial atomic solver and its
`mandacaru.pseudopotentials.oncv` module; nothing here is copied
from another pseudopotential code.

They live in this repository, not in Mandacaru itself, because of their size:
about 110 MB for the 92 files, against the 100 MB limit of a PyPI release.
The Troullier–Martins (NCPP) library, 11 MB, still ships inside the package.

## Using the datasets

Point Mandacaru at a checkout of this repository once; it creates a symbolic
link `library/oncvpsp` inside the installed package, and the loaders take it
from there:

```bash
git clone git@github.com:seixas-research/mandacaru-oncvpsp.git
python -m mandacaru.pseudopotentials.link_library --oncvpsp mandacaru-oncvpsp
```

Use `--files` to link each dataset individually instead of the directory,
`--force` to replace an existing link, `--status` to see what each family
folder serves. Alternatively set `MANDACARU_PSEUDO_PATH` to a directory that
contains this checkout as its `oncvpsp/` subfolder.

Then, in a calculation:

```python
from ase.build import molecule
from mandacaru import Mandacaru

atoms = molecule("H2O")
atoms.center(vacuum=4.0)          # the cell is the real-space box
atoms.calc = Mandacaru(method="adapt-vqe",
                       basis="ONCVPSP",
                       h=0.25)
atoms.get_total_energy()          # eV, valence-only Hamiltonian
```

The family is selected **as a basis**: `basis="ONCVPSP"` (alias `"ONCV"`), or
`basis={"name": "ONCVPSP", "size": "DZP"}` for a larger valence basis, exactly
like an all-electron family. Its guide is the *Pseudopotentials* page of the
Mandacaru manual (`docs/source/guide/pseudopotentials.md`).

## What is in a file

Each `<Symbol>.parquet` is a self-describing Mandacaru pseudopotential record
(format `mandacaru-pseudopotential`, version 2, `family = "oncvpsp"`), readable
with `mandacaru.pseudopotentials.oncv.get_oncv(symbol)` or the
generic `io.load_pseudopotential(path)`. The table holds the radial grid
(3000 points, 0.01 bohr spacing) and, per angular momentum `l`:

- two pseudo partial waves `pseudo_wave_l{l}_{0,1}` at the reference energies
  ε₁ (the bound valence eigenvalue) and ε₂ = ε₁ + 1 Ha,
- the two projectors `projector_l{l}_{0,1}` and their 2×2 Vanderbilt
  coupling matrix,
- the screened and ionic channel potentials.

The record also carries the local potential (screened and unscreened), the
pseudo valence density, the cutoff radii, the residual kinetic energies of
the optimization, the generalized-norm-conservation residuals and the
reference energies. All quantities are in atomic units (bohr, hartree);
Mandacaru converts to eV and Å at its user-facing layer.

## Construction, in brief

Following D. R. Hamann, *Phys. Rev. B* **88**, 085117 (2013):

1. all-electron LDA atom; valence partial waves at two energies per `l`;
2. pseudo partial waves as spherical-Bessel expansions inside `r_c`, matched
   in value and derivatives at `r_c`, with generalized norm conservation
   ⟨φ̃ᵢ|φ̃ⱼ⟩ = ⟨φᵢ|φⱼ⟩ inside the sphere, and the residual kinetic energy
   beyond a cutoff wavevector minimized;
3. a polynomial local potential inside `0.9 min r_c`, unscreened with the
   Hartree and LDA exchange–correlation potentials of the pseudo valence
   density;
4. projectors χᵢ = (εᵢ − T − V_loc) φ̃ᵢ with the Vanderbilt coupling
   D = B⁻¹, B_ij = ⟨φ̃ᵢ|χⱼ⟩.

Every dataset was checked on generation: the bound eigenvalue is reproduced
to better than 1e-6 Ha with no ghost state below it, the all-electron tail
is matched, the norm-conservation matrix is satisfied to 1e-8, and the
logarithmic derivatives agree with the all-electron atom at both reference
energies. Not included: nonlinear core correction, relativistic terms,
projectors above the valence `l`; the exchange–correlation functional is LDA.

## Regenerating

From the Mandacaru repository, with the `mandacaru` environment:

```python
from mandacaru.pseudopotentials.oncv import build_oncv_library
from mandacaru.pseudopotentials.io import library_elements

build_oncv_library(library_elements(92), directory="path/to/mandacaru-oncvpsp")
```

Generation takes 5 s for light elements and up to 150 s for the heaviest, about
90 minutes for the full set, with under 1 GB of memory. The 2026-09-14 set
was generated with zero failures.

## License

MIT, see `LICENSE`.
