# Simulated-Annealing-Based Back-Propagation ANN for Potential Energy Surfaces

A Fortran implementation associated with our work on using artificial neural networks to represent interatomic and intermolecular potential energy surfaces (PESs).

> **Research software status:** this repository contains the original `back-propg4.f` source supplied for this project. It is an archival research example, not a verified, maintained release. Please validate the code, compiler behavior, and numerical results before using it in research.

## About the method

The paper presents a supervised ANN approach for learning a potential-energy mapping from structure–energy data, including a simulated-annealing-based training procedure. It discusses applications to a model matrix eigenvalue/eigenvector problem, the OH diatomic potential-energy curve, and the water–formaldehyde dimer intermolecular PES.

The attached Fortran example defines a feed-forward network with **9 input values, 4 hidden units, and 4 output values**. Its input data are four hard-coded 3×3 matrices flattened into nine values, with four target values per pattern. The source perturbs one weight at a time and uses an annealing-style acceptance expression. It does not calculate analytic back-propagation gradients in the modern gradient-descent sense; the code and method description should therefore be read in the context of the cited paper.

## Files

- `back-propg4.f` — original fixed-form Fortran source, including the `ran3` random-number routine.
- `index.html` — static project landing page for GitHub Pages.
- `CITATION.cff` — citation metadata for the associated paper.

## Build and run

A Fortran compiler such as GNU Fortran (`gfortran`) is needed. From the repository root, a possible starting command is:

```bash
gfortran -std=legacy -ffixed-form -ffixed-line-length-none back-propg4.f -o backprop
./backprop
```

The command is provided as a starting point and has **not** been validated against this archival source. The program uses hard-coded patterns and declares several output files; inspect the source and confirm its behavior with your compiler before interpreting any generated output.

## Important source notes

The source is preserved as supplied rather than silently rewritten. A review before scientific reuse or public release should at least address:

- `aaexp` may be tested on an equal-error path where the current iteration has not assigned it.
- The annealing expression uses `i` in `tempt=temp0-5.0d0*(i-1)`; confirm whether the intended loop variable is `itemp`.
- Output units 4 and 74 are written without matching explicit `OPEN` statements.
- `init_err` initializes a separate random weight set from the main program's weight set.
- The `pattern` routine comments label pattern 4 as unknown, but both the pattern routine and main training loop supply/process its target values.
- The `ran3` routine is attributed in the source to *Numerical Recipes*. Check the applicable terms before redistributing that routine. No repository license is included; add an appropriate license only after confirming rights for all included code.

These are observations from the supplied file, not a claim about the intended published implementation. Compare this source against the version and procedures used for the published results.

## Publication

Mayank Dixit, Joseph Daniel, and Shankar Prasad Bhattacharyya, “Exploring artificial neural networks to model interatomic and intermolecular potential energy surfaces,” *Journal of the Indian Chemical Society* **98**(9), 100114 (2021). [https://doi.org/10.1016/j.jics.2021.100114](https://doi.org/10.1016/j.jics.2021.100114)

## GitHub Pages

This repository includes a root-level `index.html`. To publish it with GitHub Pages, push the repository to GitHub, then choose **Settings → Pages → Deploy from a branch**, select the `main` branch and `/ (root)`, and save. GitHub will provide the published URL after the deployment completes.
