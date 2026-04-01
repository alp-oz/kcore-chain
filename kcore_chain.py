"""
diamond.py

Markov chain on k-bounded partitions (R_k), defined in kcore.tex Section 3.

States: R_k = {k-bounded partitions lambda : l_i(lambda) <= k-i for all i},
        |R_k| = k!, indexed via the factorial code.

Transition probabilities (eq. 3.3 of kcore.tex):
    P(lambda, mu) = d_mu^(k) / ((|lambda|+1) * d_lambda^(k))
where d_lambda^(k) counts standard strong k-tableaux of shape c(lambda).

Dimension formula: inclusion-exclusion over raising operators R_{ij},
ported directly from dimens.m / dim.m (MATLAB). See kschur_dim().

All transition probabilities are exact fractions (fractions.Fraction).
"""

from __future__ import annotations
from math import factorial
from fractions import Fraction
from collections import defaultdict
import csv


# ---------------------------------------------------------------------------
# Partition utilities
# Partitions: tuples of positive integers in weakly decreasing order.
# Empty partition = ().
# ---------------------------------------------------------------------------

def _normalize(p: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((x for x in p if x > 0), reverse=True))


def conjugate(lam: tuple[int, ...]) -> tuple[int, ...]:
    """Transpose (conjugate) of a partition."""
    if not lam:
        return ()
    return tuple(
        sum(1 for r in range(len(lam)) if lam[r] >= c + 1)
        for c in range(lam[0])
    )


# ---------------------------------------------------------------------------
# Bijection: k-bounded partitions  <->  (k+1)-core partitions
# Maps p (core -> bounded) and c (bounded -> core), Section 2 of kcore.tex.
# ---------------------------------------------------------------------------

def bounded_to_core(lam: tuple[int, ...], k: int) -> tuple[int, ...]:
    """
    Map a k-bounded partition to its (k+1)-core.  Port of corep.m.
    Operates on the reversed (weakly-increasing) partition with a running
    shift, using the conjugate to track column heights.
    """
    if not lam:
        return ()

    h = len(lam)
    R = list(reversed(lam))        # weakly increasing (smallest part first)
    T = list(conjugate(lam))       # column lengths

    C = R[:]
    shift = 0
    count = 0                      # 0-based index into T
    s = T[0] if T else 0

    for i in range(h):
        # Condition (1-based MATLAB): i+1 + s - h + R[i] - 1 <= k
        if (i + 1) + s - h + R[i] - 1 <= k:
            C[i] = R[i] + shift
        else:
            while (i + 1) + s - h + R[i] - 1 > k:
                shift += 1
                count += 1
                s = max(T[count], h - i) if count < len(T) else (h - i)
                for z in range(i, h):
                    C[z] = R[z] + shift
                W = list(reversed(C))
                T = list(conjugate(tuple(W)))
                s = T[count] if count < len(T) else 0
            C[i] = R[i] + shift

    return _normalize(tuple(reversed(C)))


def core_to_bounded(kappa: tuple[int, ...], k: int) -> tuple[int, ...]:
    """
    Map a (k+1)-core to its k-bounded partition.
    Row i of the result has length = #{cells in row i of kappa with hook <= k}.
    """
    if not kappa:
        return ()
    n = len(kappa)
    rows = []
    for i in range(n):
        cnt = sum(
            1 for j in range(kappa[i])
            if (kappa[i] - j - 1) + sum(1 for r in range(i+1, n) if kappa[r] > j) + 1 <= k
        )
        if cnt > 0:
            rows.append(cnt)
    return _normalize(tuple(rows))


# ---------------------------------------------------------------------------
# Factorial code bijection  R_k  <->  {0, ..., k!-1}
# Port of factary.m and invfact.m.
# ---------------------------------------------------------------------------

def partition_to_index(lam: tuple[int, ...], k: int) -> int:
    """
    Encode a k-bounded partition in R_k as an integer in {0, ..., k!-1}.
    Each part p contributes k! / (k-p+1)! to the code.
    """
    n = 0
    for p in lam:
        n += factorial(k) // factorial(k - p + 1)
    return n


def index_to_partition(m: int, k: int) -> tuple[int, ...]:
    """
    Decode an integer m in {0, ..., k!-1} to a k-bounded partition in R_k.
    Inverse of partition_to_index.
    """
    if m == 0:
        return ()
    n = m % factorial(k)
    A = [0] * (k + 1)              # A[i] = # parts of size k-i+1 (1-based)
    for i in range(1, k):
        z = factorial(k) // factorial(i + 1)
        A[i + 1] = n // z
        n = n % z
    parts = []
    for i in range(1, k + 1):
        parts.extend([k - i + 1] * A[i])
    return _normalize(tuple(parts))


# ---------------------------------------------------------------------------
# Rectangle reduction: project a k-bounded partition to R_k
# Port of reduce.m.
# ---------------------------------------------------------------------------

def reduce_to_Rk(lam: tuple[int, ...], k: int) -> tuple[int, ...]:
    """
    Remove complete k-rectangles from lam to obtain its R_k representative.
    Part of size i appears mod (k-i+1) times in the reduced partition.
    """
    if not lam:
        return ()
    counts: dict[int, int] = defaultdict(int)
    for p in lam:
        counts[p] += 1
    parts = []
    for p in range(1, k + 1):
        parts.extend([p] * (counts[p] % (k - p + 1)))
    return _normalize(tuple(parts))


# ---------------------------------------------------------------------------
# Content of a (k+1)-core: weak cover step
# Port of content.m.
# ---------------------------------------------------------------------------

def _cell_residue(i: int, j: int, k: int) -> int:
    """Content residue (j - i) mod (k+1) of cell (i, j) in English notation."""
    return (j - i) % (k + 1)


def _addable_corners(kappa: tuple[int, ...], k: int, residue: int) -> list[tuple[int, int]]:
    """
    All cells that can be added to core kappa such that the result is a valid
    Young diagram and the cell has the given residue mod (k+1).
    """
    n = len(kappa)
    cells = []
    for i in range(n):
        j = kappa[i]
        above_ok = (i == 0) or (kappa[i] < kappa[i - 1])
        below_ok = (i == n - 1) or (kappa[i + 1] <= kappa[i])
        if above_ok and below_ok and _cell_residue(i, j, k) == residue:
            cells.append((i, j))
    if (n == 0 or kappa[n - 1] >= 1) and _cell_residue(n, 0, k) == residue:
        cells.append((n, 0))
    return cells


def content_partition(kappa: tuple[int, ...], k: int, residue: int) -> tuple[int, ...] | None:
    """
    Add all addable corners of given residue to core kappa simultaneously
    (one step in the weak order), then map back to a k-bounded partition.
    Returns None if no such corners exist.
    """
    cells = _addable_corners(kappa, k, residue)
    if not cells:
        return None
    lst = list(kappa)
    for (i, j) in cells:
        while len(lst) <= i:
            lst.append(0)
        lst[i] += 1
    return core_to_bounded(_normalize(tuple(lst)), k)


# ---------------------------------------------------------------------------
# Dimension  d_lambda^(k)
# Port of dim.m + dimens.m (MATLAB).
# Inclusion-exclusion over raising operators R_{ij} acting on the partition
# vector L. Uses 1-based arrays internally to match the MATLAB source exactly.
# ---------------------------------------------------------------------------

def _multinomial(n: int, parts: list[int]) -> int:
    """n! / (parts[0]! * parts[1]! * ...).  Returns 0 if any part is negative."""
    if any(v < 0 for v in parts):
        return 0
    result = factorial(n)
    for v in parts:
        result //= factorial(v)
    return result


def kschur_dim(lam: tuple[int, ...], k: int) -> int:
    """
    Compute d_lambda^(k), the number of standard strong k-tableaux of shape
    c(lambda).  Port of dim.m + dimens.m using 1-based arrays internally.
    """
    if not lam:
        return 0

    # dim.m: build multiplicity array A (1-based, length k).
    # A[k-p+1] counts parts of size p.
    A = [0] * (k + 1)
    for p in lam:
        A[k - p + 1] += 1

    p = sum(A[1:])
    if p == 0:
        return 0

    n = sum((k - i + 1) * A[i] for i in range(1, k + 1))

    # Build partition vector L and first-index array l (both 1-based).
    L = [0] * (p + 1)
    l = [0] * (k + 1)
    l[1] = 1
    count = 1
    tr = A[count]
    for i in range(1, p + 1):
        if i <= tr:
            L[i] = k - count + 1
        else:
            count += 1
            while A[count] == 0:
                l[count] = i
                count += 1
            l[count] = i
            L[i] = k - count + 1
            tr += A[count]

    # Operator counts S[j]: number of R_{ij} pairs with displacement j.
    S = [0] * k
    for i in range(1, p + 1):
        if L[i] <= k - 1:
            lim = min(k - L[i], p - i)
            for j in range(1, lim + 1):
                S[j] += 1

    s = sum(S[1:])
    if s == 0:
        return _multinomial(n, L[1:])

    # Build operator labelling arrays t and M (1-based).
    t = [0] * (s + 1)
    M = [0] * (k + 1)
    count = 1
    tr = S[count]
    for i in range(1, s + 1):
        if i <= tr:
            t[i] = count
            M[count + 1] += 1
        else:
            M[count + 2] = M[count + 1] + 1
            count += 1
            t[i] = count
            tr += S[count]

    # Inclusion-exclusion over all 2^s subsets.
    total = 0
    for mask in range(1 << s):
        K = L[:]
        z = 1
        for j in range(1, s + 1):
            if (mask >> (j - 1)) & 1:
                z = -z
                tj = t[j]
                idx_dec = l[tj + 1] + j - M[tj] + tj - 1
                idx_inc = l[tj + 1] + j - M[tj] - 1
                if 1 <= idx_dec <= p:
                    K[idx_dec] -= 1
                if 1 <= idx_inc <= p:
                    K[idx_inc] += 1
        total += z * _multinomial(n, K[1:])

    return total


# ---------------------------------------------------------------------------
# Core validity check (used internally)
# ---------------------------------------------------------------------------

def _is_core(kappa: tuple[int, ...], k: int) -> bool:
    """True iff kappa has no cell with hook length exactly k+1."""
    n = len(kappa)
    for i in range(n):
        for j in range(kappa[i]):
            arm = kappa[i] - j - 1
            leg = sum(1 for r in range(i + 1, n) if kappa[r] > j)
            if arm + leg + 1 == k + 1:
                return False
    return True


# ---------------------------------------------------------------------------
# Build the transition matrix
# ---------------------------------------------------------------------------

def markov(k: int) -> list[list[Fraction]]:
    """
    Build the k! x k! transition matrix for the Markov chain on R_k.
    Entry [i][j] = P(lambda_i -> lambda_j) as an exact Fraction.
    States are ordered by their factorial code: state i <-> index_to_partition(i, k).

    Transition probability (kcore.tex eq. 3.3):
        P(lambda, mu) = d_mu^(k) / ((|lambda|+1) * d_lambda^(k))
    where mu is the reduce_to_Rk image of a content partition of c(lambda).
    """
    size   = factorial(k)
    states = [index_to_partition(m, k) for m in range(size)]
    cores  = [bounded_to_core(lam, k) for lam in states]

    M = [[Fraction(0)] * size for _ in range(size)]

    for i, lam in enumerate(states):
        n      = sum(lam)
        d_lam  = kschur_dim(lam, k) if lam else 1
        if d_lam == 0:
            continue
        kappa  = cores[i]
        denom  = Fraction(n + 1) * d_lam

        for res in range(k + 1):
            lam_res = content_partition(kappa, k, res)
            if lam_res is None:
                continue
            d_res  = kschur_dim(lam_res, k) if lam_res else 1
            target = reduce_to_Rk(lam_res, k)
            j      = partition_to_index(target, k)
            if 0 <= j < size:
                M[i][j] += Fraction(d_res) / denom

    return M


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def save_matrix_fractions(M: list[list[Fraction]], path: str) -> None:
    """Save transition matrix as CSV of exact fractions."""
    size = len(M)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([f"B{j+1}" for j in range(size)])
        for row in M:
            writer.writerow([str(x) if x != 0 else "0" for x in row])


def save_matrix_floats(M: list[list[Fraction]], path: str) -> None:
    """Save transition matrix as CSV of floats."""
    size = len(M)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([f"B{j+1}" for j in range(size)])
        for row in M:
            writer.writerow([float(x) for x in row])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    print(f"Computing Markov chain on R_{k} ({factorial(k)} states)...")

    M = markov(k)

    import os
    os.makedirs("results", exist_ok=True)
    frac_path  = f"results/MC_{k}_frac.csv"
    float_path = f"results/MC_{k}_float.csv"
    save_matrix_fractions(M, frac_path)
    save_matrix_floats(M, float_path)
    print(f"Saved: {frac_path}, {float_path}")

    bad = sum(1 for row in M if abs(float(sum(row)) - 1.0) > 1e-10)
    if bad == 0:
        print("All row sums = 1. ✓")
    else:
        for i, row in enumerate(M):
            s = sum(row)
            if abs(float(s) - 1.0) > 1e-10:
                print(f"  WARNING: row {i} sums to {s}")
