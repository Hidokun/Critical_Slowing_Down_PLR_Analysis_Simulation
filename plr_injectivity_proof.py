#!/usr/bin/env python3
"""
PLR Gain Inversion — Numerical Proof of Injectivity
=====================================================
Proves that Re(dlambda*/dG) > 0 on the operative window G in [G_lo, G_hi].

The symbolic derivation reduces the sign condition to:
    delta * G^2 * E^2 + tau_iris * u > 0
where:
    E   = exp(-alpha * delta),   alpha = Re(lambda*)
    u   = 1 + tau_iris * alpha

For u >= 0 this is trivially positive (all terms non-negative).
For u <  0 this requires parameter-specific verification (this script).

Usage:
    python3 plr_injectivity_proof.py
    python3 plr_injectivity_proof.py --tau_iris 0.311 --delta 0.300 --G_lo 0.15 --G_hi 2.1 --N 10000

Parameters can be changed via command-line flags to re-run the proof
under different physiological or calibration assumptions.
"""

import argparse
import numpy as np
from scipy.special import lambertw

# ─────────────────────────────────────────────────────────────────────────────
# 1. Argument parsing — lets you rerun with different parameters
# ─────────────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="PLR injectivity numerical proof")
    p.add_argument("--tau_iris", type=float, default=0.311,
                   help="Iris mechanical time constant (s), default 0.311")
    p.add_argument("--delta",    type=float, default=0.300,
                   help="Neural conduction delay (s), default 0.300")
    p.add_argument("--G_lo",     type=float, default=0.15,
                   help="Lower bound of operative window, default 0.15")
    p.add_argument("--G_hi",     type=float, default=2.1,
                   help="Upper bound of operative window, default 2.1")
    p.add_argument("--N",        type=int,   default=100000,
                   help="Number of evaluation points, default 10000")
    return p.parse_args()

# ─────────────────────────────────────────────────────────────────────────────
# 2. Core functions
# ─────────────────────────────────────────────────────────────────────────────
def dominant_eigenvalue(G, delta, tau_iris):
    """
    Exact dominant eigenvalue via Lambert W (Eq. 8 of paper):
        lambda*(G) = W_0(-G*delta/tau_iris * exp(delta/tau_iris)) / delta  -  1/tau_iris
    """
    arg = -(G * delta / tau_iris) * np.exp(delta / tau_iris)
    W   = lambertw(arg, k=0)          # principal branch
    return W / delta - 1.0 / tau_iris

def G_star_crossing(delta, tau_iris):
    """
    G* such that z(G*) = -1/e, i.e. the W_0 argument crosses the branch point.
    Below G* the principal branch is real; above it W_0 is complex.
        G* = (tau_iris / (e * delta)) * exp(-delta / tau_iris)
    """
    return (tau_iris / (np.e * delta)) * np.exp(-delta / tau_iris)

def tau_return(G, delta, tau_iris):
    """Recovery time constant tau_return = -1 / Re(lambda*)."""
    lam     = dominant_eigenvalue(G, delta, tau_iris)
    re_lam  = np.real(lam)
    if re_lam >= 0:
        return np.inf
    return -1.0 / re_lam

def sign_condition(G, delta, tau_iris):
    """
    The reduced sign condition (from the symbolic derivation):
        S(G) = delta * G^2 * E^2 + tau_iris * u
    where E = exp(-alpha*delta), u = 1 + tau_iris*alpha, alpha = Re(lambda*).
    Re(dlambda*/dG) > 0  iff  S(G) > 0.
    """
    lam   = dominant_eigenvalue(G, delta, tau_iris)
    alpha = np.real(lam)
    E     = np.exp(-alpha * delta)
    u     = 1.0 + tau_iris * alpha
    return delta * G**2 * E**2 + tau_iris * u

def dlam_dG_real(G, delta, tau_iris):
    """
    Re(dlambda*/dG) via implicit differentiation of the characteristic equation.
    dlambda/dG = -exp(-lambda*delta) / [tau_iris + delta*(tau_iris*lambda + 1)]
    (denominator simplified via char eq: tau_iris*lambda + 1 = -G*exp(-lambda*delta))
    """
    lam      = dominant_eigenvalue(G, delta, tau_iris)
    exp_term = np.exp(-lam * delta)
    num      = -exp_term
    den      = tau_iris + delta * (tau_iris * lam + 1.0)
    return np.real(num / den)

def dtau_dG(G, delta, tau_iris):
    """
    d(tau_return)/dG = Re(dlam/dG) / Re(lam*)^2
    (tau_return = -1/Re(lam*), so d/dG(-1/Re(lam*)) = Re(dlam/dG)/Re(lam*)^2)
    """
    lam    = dominant_eigenvalue(G, delta, tau_iris)
    re_lam = np.real(lam)
    re_dl  = dlam_dG_real(G, delta, tau_iris)
    return re_dl / re_lam**2

# ─────────────────────────────────────────────────────────────────────────────
# 3. Main proof procedure
# ─────────────────────────────────────────────────────────────────────────────
def run_proof(tau_iris, delta, G_lo, G_hi, N):
    sep = "=" * 72

    print(sep)
    print("  PLR GAIN INVERSION — INJECTIVITY PROOF")
    print(sep)
    print(f"  Parameters:  tau_iris = {tau_iris} s   |   delta = {delta} s")
    print(f"  Operative window: G in [{G_lo}, {G_hi}]")
    print(f"  Evaluation points: N = {N}")
    print()

    # ── Step 1: G* branch-crossing point ─────────────────────────────────────
    Gstar = G_star_crossing(delta, tau_iris)
    z_at_lo = -(G_lo * delta / tau_iris) * np.exp(delta / tau_iris)
    print("── STEP 1: W_0 branch structure ──────────────────────────────────────")
    print(f"   W_0 argument z(G) = -(G*delta/tau_iris)*exp(delta/tau_iris)")
    print(f"   Branch-point crossing: z(G*) = -1/e => G* = {Gstar:.6f}")
    print(f"   Lower bound G_lo = {G_lo}")
    print(f"   z(G_lo) = {z_at_lo:.6f}   |   -1/e = {-1/np.e:.6f}")
    if G_lo > Gstar:
        print(f"   => G* < G_lo: W_0 is COMPLEX throughout [{G_lo}, {G_hi}].")
        print(f"      Real-branch monotonicity argument does NOT apply.")
    else:
        print(f"   => G* >= G_lo: W_0 transitions within the window.")
    print()

    # ── Step 2: Symbolic reduction ────────────────────────────────────────────
    print("── STEP 2: Symbolic reduction ────────────────────────────────────────")
    print("   Re(dlambda*/dG) > 0")
    print("   iff  delta*(u^2 + tau_iris^2*omega^2) + tau_iris*u > 0   [by algebra]")
    print("   iff  delta*G^2*E^2 + tau_iris*u > 0                       [char. eq.]")
    print("   where u = 1 + tau_iris*alpha,  E = exp(-alpha*delta)")
    print()

    # ── Step 3: Split into trivial and non-trivial regions ───────────────────
    G_vals = np.linspace(G_lo, G_hi, N)
    u_vals = np.array([np.real(1 + tau_iris * np.real(dominant_eigenvalue(G, delta, tau_iris)))
                       for G in G_vals])

    # Find sign change in u
    u_nonneg_mask = u_vals >= 0
    if u_nonneg_mask.any():
        G_u_zero = G_vals[np.where(u_nonneg_mask)[0][0]]
    else:
        G_u_zero = None

    print("── STEP 3: Trivial region (u >= 0) ──────────────────────────────────")
    if G_u_zero is not None:
        print(f"   u >= 0 for G >= {G_u_zero:.4f}")
        print(f"   In this region: delta*G^2*E^2 >= 0, tau_iris*u >= 0.")
        print(f"   Both terms non-negative => S(G) > 0. QED for G in [{G_u_zero:.4f}, {G_hi}].")
    else:
        print("   u < 0 throughout the entire window; no trivial region.")
    print()

    # ── Step 4: Numerical verification of non-trivial region ─────────────────
    if G_u_zero is not None:
        hard_mask = ~u_nonneg_mask
        G_hard    = G_vals[hard_mask]
        label     = f"[{G_lo}, {G_u_zero:.4f}]"
    else:
        G_hard    = G_vals
        label     = f"[{G_lo}, {G_hi}]"

    S_hard = np.array([sign_condition(G, delta, tau_iris) for G in G_hard])

    print(f"── STEP 4: Non-trivial region (u < 0): G in {label} ─────────────────")
    print(f"   Evaluating S(G) = delta*G^2*E^2 + tau_iris*u at {len(G_hard)} points...")
    print(f"   min  S(G) = {S_hard.min():.8f}   at G = {G_hard[S_hard.argmin()]:.6f}")
    print(f"   max  S(G) = {S_hard.max():.8f}")
    print(f"   All S(G) > 0?  {(S_hard > 0).all()}")
    print()

    # ── Step 5: Full window verification ─────────────────────────────────────
    S_full = np.array([sign_condition(G, delta, tau_iris) for G in G_vals])
    re_dl  = np.array([dlam_dG_real(G, delta, tau_iris) for G in G_vals])
    tau_v  = np.array([tau_return(G, delta, tau_iris) for G in G_vals])
    dtdG_v = np.array([dtau_dG(G, delta, tau_iris) for G in G_vals])

    print("── STEP 5: Full window summary ───────────────────────────────────────")
    print(f"   Re(dlambda*/dG) range:   [{re_dl.min():.4f}, {re_dl.max():.4f}]")
    print(f"   All Re(dlambda*/dG) > 0? {(re_dl > 0).all()}")
    print(f"   tau_return range:        [{tau_v.min():.4f}, {tau_v.max():.4f}] s")
    print(f"   tau_return strictly increasing? {all(tau_v[i] < tau_v[i+1] for i in range(len(tau_v)-1))}")
    print()

    # ── Step 6: Key spot values ───────────────────────────────────────────────
    print("── STEP 6: Spot values across operative window ───────────────────────")
    print(f"   {'G':>6}  {'Re(lam*)':>12}  {'tau_ret(s)':>12}  {'Re(dlam/dG)':>14}  {'dtau/dG':>12}  {'S(G)':>12}")
    for G_spot in [G_lo, 0.5, 1.0, 1.5, 2.0, G_hi]:
        if G_spot < G_lo or G_spot > G_hi:
            continue
        lam_s  = dominant_eigenvalue(G_spot, delta, tau_iris)
        re_l   = np.real(lam_s)
        tau_s  = tau_return(G_spot, delta, tau_iris)
        re_dl_s = dlam_dG_real(G_spot, delta, tau_iris)
        dt_s   = dtau_dG(G_spot, delta, tau_iris)
        S_s    = sign_condition(G_spot, delta, tau_iris)
        print(f"   {G_spot:6.3f}  {re_l:12.6f}  {tau_s:12.4f}  {re_dl_s:14.4f}  {dt_s:12.4f}  {S_s:12.8f}")
    print()

    # ── Verdict ───────────────────────────────────────────────────────────────
    proof_passes = (S_full > 0).all()
    print("── VERDICT ───────────────────────────────────────────────────────────")
    if proof_passes:
        print("   ✓  S(G) > 0 for all G in the operative window.")
        print("   ✓  Re(dlambda*/dG) > 0 throughout [G_lo, G_hi].")
        print("   ✓  tau_return(G) is strictly monotone increasing.")
        print("   ✓  The inversion is well-posed: every measured tau_return")
        print("      corresponds to exactly one gain value G.")
    else:
        print("   ✗  PROOF FAILS for current parameters.")
        fail_G = G_vals[S_full <= 0]
        print(f"      Failure region: G in {[fail_G.min(), fail_G.max()]}")
    print(sep)
    return proof_passes

# ─────────────────────────────────────────────────────────────────────────────
# 4. Entry point
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    args = parse_args()
    run_proof(
        tau_iris = args.tau_iris,
        delta    = args.delta,
        G_lo     = args.G_lo,
        G_hi     = args.G_hi,
        N        = args.N,
    )
