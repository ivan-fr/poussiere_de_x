import matplotlib.pyplot as plt
import numpy as np
import random
import os

# T4 Complex Simulation
def _complex_pandrosion_T4(x_c: complex, p: float = 3.0) -> complex:
    if x_c == 1.0 + 0j or x_c == 0.0 + 0j: return 0j
    s_0 = 1.0 - (x_c - 1.0) / (x_c * p)
    
    def S_p(s):
        if s == 1.0 + 0j: return p
        return (1.0 - s**p) / (1.0 - s)
        
    def h(s):
        sp = S_p(s)
        if sp == 0j: return s
        return 1.0 - (x_c - 1.0) / (x_c * sp)
        
    s_1 = h(s_0)
    s_2 = h(s_1)
    
    diff_s1_s0 = s_1 - s_0
    hat_lambda = (s_2 - s_1) / diff_s1_s0 if abs(diff_s1_s0) > 1e-12 else 0j
    
    denom_t2 = (s_2 - 2*s_1 + s_0)
    T_2 = s_0 - (diff_s1_s0**2) / denom_t2 if abs(denom_t2) > 1e-12 else s_0
    
    s_3 = h(T_2)
    denom_t4 = hat_lambda - 1.0
    T_4 = T_2 - (s_3 - T_2) / denom_t4 if abs(denom_t4) > 1e-12 else T_2
    
    residual_gap = s_0 - T_4
    return residual_gap, hat_lambda

def generate_figure():
    np.random.seed(42)
    # Simulate a normalized geometric walk representing (x/A) fractions
    steps = 100
    prices = [1.0]
    for _ in range(steps):
        change = np.random.normal(0, 0.005)
        # Induce a structural divergence dip in the middle
        if 40 < _ < 60: change -= 0.01 
        prices.append(prices[-1] + change)
        
    baselines, uppers, lowers, lambdas = [], [], [], []
    
    for i in range(len(prices)):
        if i < 15:
            baselines.append(np.nan)
            uppers.append(np.nan)
            lowers.append(np.nan)
            lambdas.append(np.nan)
            continue
            
        slice_p = prices[i-15:i+1]
        p_now = slice_p[-1]
        p_past = slice_p[0]
        
        baseline = np.mean(slice_p)
        ratio = p_now / p_past if p_past > 0 else 1.0
        
        # Calculate
        gap, lam = _complex_pandrosion_T4(complex(ratio, 0))
        vol = 3.0 * np.sqrt(abs(gap))
        
        baselines.append(baseline)
        uppers.append(baseline * (1 + vol))
        lowers.append(baseline * (1 - vol))
        lambdas.append(lam.real * 100.0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]})
    
    ax1.plot(prices, label='Price Context (x)', color='white')
    ax1.plot(baselines, label='Baseline Anchor (A)', color='yellow', linestyle='--')
    ax1.plot(uppers, label=r'Upper Dust Band ($\sqrt{\epsilon}$)', color='red', alpha=0.7)
    ax1.plot(lowers, label=r'Lower Dust Limit ($\sqrt{\epsilon}$)', color='green', alpha=0.7)
    ax1.set_title("Pandrosion Mathematical Structural Bounds")
    ax1.set_facecolor('#1e1e1e')
    ax1.legend()
    ax1.grid(True, alpha=0.2)
    
    colors = ['green' if l > 0 else 'red' for l in lambdas]
    ax2.bar(range(len(lambdas)), lambdas, color=colors)
    ax2.set_facecolor('#1e1e1e')
    ax2.set_title(r"Structural Eigenvalue ($\hat{\lambda}$)")
    fig.patch.set_facecolor('#121212')
    
    for ax in [ax1, ax2]:
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.tick_params(colors='white')
        
    plt.tight_layout()
    output_path = os.path.join(os.path.dirname(__file__), '../figures/pandrosion_t4_bounds.png')
    plt.savefig(output_path, dpi=300, facecolor=fig.get_facecolor(), bbox_inches='tight')
    print(f"Generated Figure at {output_path}")

if __name__ == "__main__":
    generate_figure()
