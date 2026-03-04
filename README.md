[![DOI](https://zenodo.org/badge/1108516148.svg)](https://doi.org/10.5281/zenodo.17807157)
[![Hippocratic License HL3-BDS-CL-FFD-LAW-MEDIA-MIL-SV](https://img.shields.io/static/v1?label=Hippocratic%20License&message=HL3-BDS-CL-FFD-LAW-MEDIA-MIL-SV&labelColor=5e2751&color=bc8c3d)](https://firstdonoharm.dev/version/3/0/bds-cl-ffd-law-media-mil-sv.html)

# DiscConv
Newtonian potential of the indicator function of a disc intersection $● := \mathrm{B}[0;1]\cap \mathrm{B}[x;\varepsilon]$

### How do I run tests? ###
 
1. Create and activate a Python environment (e.g. via `venv` or `conda`) with Python ≥ 3.10.

2. Install the required Python libraries:
   ```bash
   pip install numpy scipy sympy mpmath matplotlib
   ```

The experiments in the note were run using:

* `Python 3.10.14`
* `NumPy 2.2.5`
*
* `SymPy 1.3.3`
* `mpmath 1.3.0`

3. Open the Jupyter notebooks in numerical order (`1 - ...` to `8 - ...`) and execute all cells.
   Each notebook reproduces and validates a specific part of the analyses.

---

## Notebooks in this folder

The numerical experiments in the note are organised into eight notebooks:

1. [`1 - Disc Integral.ipynb`](1%20-%20Disc%20Integral.ipynb)
   Derives and tests the closed-form expression
   $K \ast \mathbf{1}_{●}(x) = \frac{1}{4}\big(|x|^2 - 1\big)$
   for the Newtonian potential of the unit disc.

   * Provides an alternative derivation of equation (5) using a log–cosine integral.
   * Compares the analytic formula against direct SciPy quadrature.

2. [`2 - Convolution inside disc.ipynb`](2%20-%20Convolution%20inside%20disc.ipynb)
   Focuses on the *nested–intersection* regime $a \in [0,1-\varepsilon]$ and on the branch where the intersection is itself a disc.

   * Numerically checks equation (7).
   * Symbolically derives the link between the dilogarithm and the log–cosine integrand using SymPy.
   * Numerically validates the intermediate identities (e.g. (17) and (18)) against quadrature.
   * Implements the representation (13a) and compares it with several quadrature methods.

3. [`3 - Direct and Error.ipynb`](3%20-%20Direct%20and%20Error.ipynb)
   Studies the *partial–overlap* regime $|x|\in (1-\varepsilon,1+\varepsilon]$ and the alternative expression for $F_\varepsilon$.

   * Symbolically verifies identities such as (25) and (26).
   * Implements the closed form (21) and tests it against:

     * the original integral definition,
     * the representation (13a), and
     * different quadrature schemes.
   * Visualises scaling and numerical errors across the overlap region.

4. [`4 - Tests H - Branch 1.ipynb`](4%20-%20Tests%20H%20-%20Branch%201.ipynb)
   High–precision tests for the reparametrised quantity $H(\lambda,\varepsilon)$ on the branch
   $a^2 \leq 1 + \varepsilon^2$,
   corresponding to one geometric regime in Lemma 4.

   * Evaluates (H) and its asymptotic expansion in extended precision.
   * Produces the error curves appearing in Figure 3 for this branch.

5. [`5 - Tests H - Branch 2.ipynb`](5%20-%20Tests%20H%20-%20Branch%202.ipynb)
   Companion to notebook 4 for the second branch $a^2 > 1 + \varepsilon^2$.

   * Evaluates $H$ and its asymptotic expansion in extended precision on this branch.
   * Completes the error plots of Figure 3 and confirms that the maximum error behaves like $\varepsilon \times 10^{-2}$.

6. [`6 - H - High double.ipynb`](6%20-%20H%20-%20High%20double.ipynb)
   Double–precision implementation of $H(\lambda,\varepsilon)$ based on the asymptotic expansion in Lemma 4.

7. [`7 - G - Computation.ipynb`](7%20-%20G%20-%20Computation.ipynb)
   Computes $G$ in extended precision for the particular pairs $(a,\Phi(\varphi_\varepsilon(a)))$, $(a,\Phi(\alpha))$, and $(a,\Phi(\varphi_\varepsilon(a)+pi))$ focusing on their relevant branches. Particular care is given on the special point $a = 1$, where the two branches of $F_\varepsilon$ meet. Additionally, all tests are validated using extended-precision quadrature.

9. [`8 - E - High double.ipynb`](8%20-%20E%20-%20High%20double.ipynb)
   Double–precision implementation of the full convolution $E_\varepsilon$ using the asymptotic information on $F_\varepsilon$.

   * For $\varepsilon \leq 10^{-5}$, replaces direct evaluation of $F_\varepsilon$ by the leading asymptotic contribution given by $h_1^{(1)}$.
   * Reproduces the high–precision plots of Figure 5 (scaled by $\varepsilon^2 \log \varepsilon^2$) down to $\varepsilon = 10^{-30}$.
   * Compares extended–precision and double–precision versions of the scaled quantity, showing agreement to about $10^{-7}$ for $\varepsilon \gtrsim 10^{-10}$ and that the apparent loss of digits in the scaled plots is due solely to ill-conditioning of the scaling factor.

---

### Additional code

Beyond the notebooks, the repository also contains:

* [`Conv_Intersection.py`](Conv_Intersection.py) A **Python module** for evaluating $E_\varepsilon$ as a radial function, intended for use in nonlocal PDE solvers and error–analysis experiments.
* [`Conv_Disc_Intersection.m`](Conv_Intersection.m) A **Matlab implementation** of the same radial evaluation, mirroring the Python interface for ease of integration into existing Matlab workflows.

These radial routines are the ones referred to at the end of the note and are suitable for direct use in applications once the parameters $\varepsilon$ and the radius $a = |x|$ are specified.

To run the latter, you can test
```
% Restricted-range example for E_np
% a ∈ [1-2ε, 1+2ε], ε = 1e-2

epsilon = 1e-2;

% Dense sampling near the transition
a = linspace(1 - 2*epsilon, 1 + 2*epsilon, 1e4);

% Exact version
E_exact = Conv_Disc_Intersection(a, epsilon, false);

% Asymptotic version (optional, to compare)
E_asymp = Conv_Disc_Intersection(a, epsilon, true);

% Plot
figure;
plot(a, E_exact, 'LineWidth', 1.5); hold on;
plot(a, E_asymp, '--', 'LineWidth', 1.2);
xlabel('a');
ylabel(sprintf('E(a; \\epsilon = %.0e)', epsilon));
title('E\_np near a \approx 1 \pm \epsilon');
legend('Exact', 'Asymptotic', 'Location', 'best');
grid on;
```



---

### Reference

If you use this code in your work, please cite:

> A. Miniguano-Trujillo, *The Newtonian kernel at the intersection of two discs*, 2025.
