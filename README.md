# kcore-chain

Computes the Markov chain on **k-bounded partitions** (the set R_k) defined in:

> The k-Plancherel measure and a Finite Markov Chain (with Svante Linusson), https://arxiv.org/abs/2512.24346

## Background

States are k-bounded partitions in R_k (|R_k| = k!), indexed via the factorial code.
Transition probabilities follow the equation in Sect. 3.3 of the paper:

```
P(λ, μ) = d_μ^(k) / ((|λ|+1) · d_λ^(k))
```

where d_λ^(k) counts standard strong k-tableaux of shape c(λ), and μ is obtained
from λ via a weak cover step on its (k+1)-core followed by rectangle reduction.
All probabilities are computed as exact fractions.

## Usage

```bash
python kcore_chain.py <k>
```

Default is k=3. Output is saved to `results/`:

- `results/MC_<k>_frac.csv` — transition matrix as exact fractions
- `results/MC_<k>_float.csv` — transition matrix as floats

## Data

Exact values obtained via separate computation for k = 3, 4, 5, 6 are available in the `Data/` folder:

- `Data/transition_matrices/transition_matrix_k<k>.csv` — transition matrix M_k as exact fractions
- `Data/stationary_distributions/stationary_distribution_k<k>.csv` — stationary distribution π_k as exact fractions

## Dependencies

Python standard library only (`fractions`, `math`, `collections`, `csv`).
