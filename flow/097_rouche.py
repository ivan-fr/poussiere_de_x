"""
PAPER: 097 (canonical: 87pandrosion_rouche.pdf)
TITLE: Strong Rouché via Pandrosion Winding
STATUS: proved (Rouché 1862, Pandrosion-form)
DEPENDS: 016, 086

THEORY
======

ROUCHÉ: If |P(z) - Q(z)| < |Q(z)| on contour gamma, then P and Q have the
same number of roots inside gamma.

PANDROSION FORM: ind(P, gamma) = ind(Q, gamma) when the perturbation is
small (does not cross any zero of Q). This is the WINDING-NUMBER stability
under perturbation.

VERIFICATION
============

  1. Rouché on test pairs.
  2. Index preservation under epsilon perturbation.
"""
from __future__ import annotations
import numpy as np


def winding(P, center, R, n_pts=512):
    thetas = 2 * np.pi * np.arange(n_pts) / n_pts
    z = center + R * np.exp(1j * thetas)
    Pp = np.polyder(P)
    F = np.polyval(Pp, z) / np.polyval(P, z)
    dz = 1j * R * np.exp(1j * thetas)
    integral = np.mean(F * dz) * 2 * np.pi
    return integral / (2 * np.pi * 1j)


def main():
    print("=" * 80)
    print("PAPER 97 — Strong Rouché via Pandrosion winding")
    print("=" * 80)

    print("\n[1] Q = z^4 + 5: 4 roots, all far from 0")
    Q = np.array([1.0, 0, 0, 0, 5])
    R = 2.0
    w_Q = winding(Q, 0, R)
    print(f"  ind(Q, |z|=2) = {w_Q.real:.4f}")

    print("\n[2] P = Q + epsilon * z: small perturbation, same index")
    for eps in [0.01, 0.5, 1.0, 5.0]:
        P = np.array([1.0, 0, 0, eps, 5])
        # Check Rouché condition: |P - Q| = eps |z| <= 2 eps on |z|=2; |Q| >= 5 - 2^4 = -11... need more care
        # Actually: |Q(z)| at |z|=2 lower bound; check empirically
        thetas = np.linspace(0, 2*np.pi, 200)
        z = 2 * np.exp(1j * thetas)
        Q_vals = np.abs(np.polyval(Q, z))
        diff_vals = np.abs(np.polyval(P - Q, z))
        rouche_ok = (diff_vals < Q_vals).all()
        w_P = winding(P, 0, R)
        print(f"  eps = {eps}: Rouché holds = {rouche_ok}, ind(P) = {w_P.real:.4f}")

    print("\n[3] Multiplicity preservation")
    P = np.array([1.0, 0, 0, 0, -1])  # z^4 - 1
    for R in [0.5, 1.5, 2.0]:
        w = winding(P, 0, R)
        n_inside = sum(1 for r in np.roots(P) if abs(r) < R)
        print(f"  R = {R}: ind = {w.real:.4f}, roots inside = {n_inside}")


if __name__ == "__main__":
    main()
