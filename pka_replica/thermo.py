"""Thermodynamics: microstate free energies -> partition function -> observables.

Implements eqs. (1)-(4) of the Uni-pKa / Starling formalism with dimensionless
free energies g = beta*G (energies from the backend are in pKa*log(10) units,
i.e. already divided by RT):

    g_i(pH)      = g_i + z_i * ln(10) * pH                      (eq 3)
    w_i(pH)      = exp(-g_i(pH)) / Z(pH)                        (eq 4)
    G_deg(z)     = -ln sum_{i in z} exp(-g_i)                   (degenerate macrostate)
    pKa(z->z-1)  = (G_deg(z) - G_deg(z-1)) / ln(10)             (eq 2, adjacent charges)

    <z>(pH)      = sum_i w_i(pH) * z_i
    pI           = pH where <z>(pH) = 0
    logD(pH)     = sum_i w_i(pH) * logP_i   (population-weighted Crippen logP,
                                             exactly as in Rowan's macropka docs)
"""

from __future__ import annotations

import math

LN10 = math.log(10.0)


def ph_dependent(g: float, charge: int, pH: float) -> float:
    """Dimensionless pH-dependent free energy (eq 3)."""
    return g + charge * LN10 * pH


def microstate_weights(
    energies: list[float], charges: list[int], pH: float
) -> list[float]:
    """Boltzmann populations of every microstate at a given pH (eq 4)."""
    gs = [ph_dependent(g, z, pH) for g, z in zip(energies, charges)]
    mn = min(gs)
    exps = [math.exp(-(g - mn)) for g in gs]
    Z = sum(exps)
    return [e / Z for e in exps]


def macrostate_free_energy(energies: list[float], charges: list[int], z: int) -> float | None:
    """G_deg(z) = -ln sum exp(-g_i) over microstates carrying charge z (log-sum-exp)."""
    gs = [g for g, zi in zip(energies, charges) if zi == z]
    if not gs:
        return None
    mn = min(gs)
    return mn - math.log(sum(math.exp(-(g - mn)) for g in gs))


def macro_pka_values(
    energies: list[float], charges: list[int]
) -> list[dict]:
    """Adjacent-charge macro-pKas (eq 2), highest charge first.

    The z_hi -> z_lo transition is a deprotonation (loss of H+), so
    pKa = (G_deg(z_lo) - G_deg(z_hi)) / ln(10).
    """
    zs = sorted({z for z in charges}, reverse=True)
    out: list[dict] = []
    for z_hi, z_lo in zip(zs, zs[1:]):
        if z_hi - z_lo != 1:
            continue
        g_hi = macrostate_free_energy(energies, charges, z_hi)
        g_lo = macrostate_free_energy(energies, charges, z_lo)
        assert g_hi is not None and g_lo is not None
        out.append(
            {
                "initial_charge": z_hi,
                "final_charge": z_lo,
                "pKa": round((g_lo - g_hi) / LN10, 3),
            }
        )
    return out


def average_charge(
    energies: list[float], charges: list[int], pH: float
) -> float:
    w = microstate_weights(energies, charges, pH)
    return sum(wi * zi for wi, zi in zip(w, charges))


def isoelectric_point(
    energies: list[float], charges: list[int],
    min_pH: float = 0.0, max_pH: float = 14.0,
) -> float | None:
    """pH where <z> crosses zero, refined by bisection on a 0.1 grid."""
    grid = [min_pH + 0.1 * i for i in range(int(round((max_pH - min_pH) / 0.1)) + 1)]
    prev = None
    bracket = None
    for pH in grid:
        z = average_charge(energies, charges, pH)
        if prev is not None and (prev[1] > 0 >= z or prev[1] < 0 <= z):
            bracket = (prev[0], pH)
            break
        prev = (pH, z)
    if bracket is None:
        return None
    lo, hi = bracket
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if average_charge(energies, charges, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round(0.5 * (lo + hi), 3)


def population_profile(
    energies: list[float], charges: list[int],
    min_pH: float = 0.0, max_pH: float = 14.0, step: float = 0.1,
) -> list[list]:
    """[[pH, [w_0, w_1, ...]], ...] — Rowan's microstate_weights_by_pH (upper bound excluded)."""
    n = int(round((max_pH - min_pH) / step))
    out = []
    for i in range(n):  # upper endpoint excluded, as in Rowan
        pH = round(min_pH + step * i, 3)
        w = microstate_weights(energies, charges, pH)
        out.append([pH, [round(x, 6) for x in w]])
    return out


def logd_profile(
    energies: list[float], charges: list[int], logps: list[float],
    min_pH: float = 0.0, max_pH: float = 14.0, step: float = 0.1,
) -> list[list]:
    """logD(pH) = sum_i w_i(pH) * logP_i (Crippen logP per microstate)."""
    n = int(round((max_pH - min_pH) / step))
    out = []
    for i in range(n):
        pH = round(min_pH + step * i, 3)
        w = microstate_weights(energies, charges, pH)
        out.append([pH, round(sum(wi * pi for wi, pi in zip(w, logps)), 3)])
    return out
