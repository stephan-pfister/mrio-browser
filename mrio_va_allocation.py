"""
MRIO Value Added Allocation
============================
For each sector, output is split into:
  - VA share (f = va / x):         satellite allocated directly to the producing sector
  - Intermediate share (1 - f):    satellite traced through the full Leontief inverse

Total satellite impact = direct + indirect
"""

import numpy as np


def compute_leontief(Z: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Compute Leontief inverse L = (I - A)^{-1}."""
    n = Z.shape[0]
    A = Z / x[np.newaxis, :]        # technical coefficients (n, n)
    L = np.linalg.inv(np.eye(n) - A)
    return L


def va_allocation(
    Z: np.ndarray,
    Y: np.ndarray,
    S: np.ndarray,
    va: np.ndarray,
) -> dict:
    """
    Allocate satellite accounts to final demand categories using VA shares.

    Parameters
    ----------
    Z  : (n, n)  Intermediate transaction matrix
    Y  : (n, m)  Final demand matrix
    S  : (k, n)  Satellite matrix  (e.g. emissions, employment)
    va : (n,)    Value added vector

    Returns
    -------
    dict with keys:
        'direct'   (k, m)  Satellite from VA-share, allocated directly
        'indirect' (k, m)  Satellite from intermediate-input share, via Leontief
        'total'    (k, m)  Sum of direct + indirect
        'f'        (n,)    VA share per sector
        'L'        (n, n)  Leontief inverse
        'e'        (k, n)  Satellite intensity (per unit output)
    """
    n = Z.shape[0]

    # Total output: row sums of Z plus row sums of Y
    x = Z.sum(axis=1) + Y.sum(axis=1)          # (n,)

    # VA share per sector
    f = va / x                                   # (n,)

    # Satellite intensity
    e = S / x[np.newaxis, :]                    # (k, n)

    # Leontief inverse
    L = compute_leontief(Z, x)                  # (n, n)

    # Direct: VA-share portion allocated without supply chain tracing
    direct = (e * f[np.newaxis, :]) @ Y         # (k, m)

    # Indirect: intermediate-input share traced through full supply chain
    indirect = (e * (1 - f)[np.newaxis, :]) @ (L @ Y)   # (k, m)

    total = direct + indirect                   # (k, m)

    return {
        "direct": direct,
        "indirect": indirect,
        "total": total,
        "f": f,
        "L": L,
        "e": e,
        "x": x,
    }


# ---------------------------------------------------------------------------
# Example usage
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    rng = np.random.default_rng(42)

    # Small 3-region × 2-sector MRIO (n = 6 sectors total)
    n = 6
    m = 3   # final demand categories
    k = 2   # satellite accounts (e.g. CO2, employment)

    # Random intermediate matrix (keep values small so A is well-behaved)
    Z = rng.uniform(10, 50, (n, n))
    Y = rng.uniform(20, 100, (n, m))

    # Total output
    x = Z.sum(axis=1) + Y.sum(axis=1)

    # Value added: remainder after subtracting intermediate inputs (column sums of Z)
    va = x - Z.sum(axis=0)
    va = np.maximum(va, 0)   # ensure non-negative

    # Satellite matrix (e.g. tons CO2, number of workers)
    S = rng.uniform(5, 30, (k, n))

    result = va_allocation(Z, Y, S, va)

    print("VA shares per sector (f):")
    print(np.round(result["f"], 3))

    print("\nDirect satellite by final demand category:")
    print(np.round(result["direct"], 2))

    print("\nIndirect satellite by final demand category:")
    print(np.round(result["indirect"], 2))

    print("\nTotal satellite by final demand category:")
    print(np.round(result["total"], 2))

    # Sanity: compare with standard Leontief result  e @ L @ Y
    e = result["e"]
    standard = e @ result["L"] @ Y
    print("\nStandard Leontief result (e @ L @ Y) for comparison:")
    print(np.round(standard, 2))
