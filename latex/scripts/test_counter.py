import numpy as np
import math

def eval_poly(roots, z):
    return np.prod(z - roots)

def Q_diff(roots, z0, z):
    if abs(z - z0) < 1e-18:
        mat = z0 - roots
        prod_tot = np.prod(mat)
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.sum(np.where(mat != 0, prod_tot / mat, 0))
    return (eval_poly(roots, z) - eval_poly(roots, z0)) / (z - z0)

def sigma_log_and_Q(roots, z0, z):
    d = len(roots)
    dz = z - roots
    dz0 = z0 - roots
    Pz = np.prod(dz)
    Pz0 = np.prod(dz0)
    
    if np.any(np.abs(dz) < 1e-18) or np.any(np.abs(dz0) < 1e-18):
        Rk_z = np.array([np.prod(np.delete(dz, k)) for k in range(d)])
        Rk_z0 = np.array([np.prod(np.delete(dz0, k)) for k in range(d)])
    else:
        Rk_z = Pz / dz
        Rk_z0 = Pz0 / dz0
        
    if abs(z - z0) < 1e-18:
        Qks = np.zeros(d, dtype=complex)
        for k in range(d):
            for m in range(d):
                if m == k: continue
                t = 1.0 + 0j
                for j in range(d):
                    if j == k or j == m: continue
                    t *= (z0 - roots[j])
                Qks[k] += t
    else:
        Qks = (Rk_z - Rk_z0) / (z - z0)
        
    Q = (Pz - Pz0) / (z - z0) if abs(z - z0) > 1e-18 else Q_diff(roots, z0, z)
    if abs(Q) < 1e-80:
        return None, Q
    absQks = np.abs(Qks)
    absQ = abs(Q)
    mask = absQks > 1e-80
    if not mask.any():
        return None, Q
    logs = np.log(absQks[mask] / absQ)
    return float(np.sum(logs)), Q

def pandrosion_T3_step(roots, z0, z):
    P0 = eval_poly(roots, z0)
    Q = Q_diff(roots, z0, z)
    if abs(Q) < 1e-60: return None
    s1 = z0 - P0 / Q
    Q2 = Q_diff(roots, z0, s1)
    if abs(Q2) < 1e-60: return s1
    s2 = z0 - P0 / Q2
    den = s2 - 2 * s1 + z
    if abs(den) < 1e-60: return s1
    return z - (s1 - z) ** 2 / den

roots = np.array([-0.96731984-0.06459744j,  0.22987706-0.45925592j,  0.06489314-0.22842954j])
z0 = -0.6496697856287815-0.28032195127130505j

print(f"z0: {z0}")
z_c = z0
print(f"P(z0) = {eval_poly(roots, z0)}")
for i in range(1):
    z_new = pandrosion_T3_step(roots, z0, z_c)
    z_c = z_new
    print(f"z_{i+1}: {z_c}")
    print(f"P(z_{i+1}) = {eval_poly(roots, z_c)}")

sigma, _ = sigma_log_and_Q(roots, z0, z_c)
print(f"Sigma_log: {sigma}")
print(f"log|P_K/P_0|: {math.log(abs(eval_poly(roots, z_c)) / abs(eval_poly(roots, z0)))}")
