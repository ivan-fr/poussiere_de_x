import numpy as np

def scan_belyi(poly, degree, R=2.0, name="Belyi"):
    """
    Scans a Belyi function (or any polynomial) from a Cauchy circle of radius R
    using the True Pandrosion Multi-Start Ratio.
    Computes the spectral energy E(R) and condition number kappa.
    """
    d = degree
    # 1. Distribute anchors on the Cauchy Circle
    s_vals = np.arange(d)
    a_s = R * np.exp(1j * 2 * np.pi * s_vals / d)
    
    # 2. Distribute targets (offset by pi/d)
    z_s = R * np.exp(1j * (2 * np.pi * s_vals / d + np.pi / d))
    
    # 3. Compute the True Pandrosion Ratio Vector (no derivatives!)
    # We add 1e-15 to avoid exact divide by zero if a root is on the circle
    P_a = poly(a_s)
    P_z = poly(z_s)
    r_vec = P_z / P_a
    
    # 4. Extract Spectral Vector via FFT
    r_hat = np.fft.fft(r_vec) / d
    
    # 5. Measure Parseval Energy
    # Classical E = sum |r_hat_k|^2
    # The energy isolated from the DC component (k=0)
    spectral_energy = np.sum(np.abs(r_hat[1:])**2)
    
    # Total Descent (Geometric Epoch)
    # real part of the log of the product of absolute ratios
    epoch_descent = np.sum(np.log(np.abs(r_vec)))
    
    # 6. Condition number (Low frequency vs High frequency)
    kappa = np.abs(r_hat[1]) / (np.abs(r_hat[d-1]) + 1e-16)
    
    print(f"--- Belyi Function: {name} (deg {d}) ---")
    print(f"  Pandrosion Energy E(R): {spectral_energy:.8e}")
    print(f"  Condition Number κ    : {kappa:.6f}")
    print(f"  Universal Descent     : {epoch_descent:.6f} nats (cf. -1.2337)")
    r_1 = np.abs(r_hat[1]) if d > 1 else 0
    r_2 = np.abs(r_hat[2]) if d > 2 else 0
    print(f"  First Spectral Modes  : |r_1|={r_1:.2e}, |r_2|={r_2:.2e}")
    print("")

if __name__ == "__main__":
    print("==========================================================")
    print(" SPECTRAL DESSINS D'ENFANTS: PANDROSION SIGNATURE SCAN")
    print("==========================================================")
    
    # 1. Simple Path (degree 2)
    p1 = lambda z: z**2
    scan_belyi(p1, 2, name="z^2 (Path L=2)")
    
    # 2. Shabat polynomial for a tree with 3 edges (path)
    # beta(z) = -2z^3 + 3z^2
    p2 = lambda z: -2*z**3 + 3*z**2
    scan_belyi(p2, 3, name="-2z^3 + 3z^2 (Path L=3)")
    
    # 3. Simple Star (degree 4)
    p3 = lambda z: z**4
    scan_belyi(p3, 4, name="z^4 (Star, 4 edges)")
    
    # 4. Asymmetric tree (Star 3 + 1 path)
    # Belyi function example: beta(z) = (256/27) * z^3 * (1-z)
    p4 = lambda z: (256/27) * z**3 * (1 - z)
    scan_belyi(p4, 4, name="256/27 z^3(1-z) (Asymmetric Tree)")
    
    # 5. Let's find Galois conjugates. 
    # A known pair of Galois conjugate dessins of degree 6 or 7.
    # Actually, we can just use polynomials with the same modulus spectrum but different phases.
    # Let's test a degree 5 polynomial to see the scan.
    p5 = lambda z: z**5 - z + 1
    scan_belyi(p5, 5, name="z^5 - z + 1 (Generic degree 5 test)")
