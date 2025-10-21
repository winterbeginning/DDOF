import re
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

THERMO_FILE = "thermoPolynomial"
MECH = "chem.yaml"
SPECIES = ["O2", "N2"]
SPLIT_T = 1000.0

def parse_paren_nums(s):
    m = re.search(r'\(([^)]*)\)', s)
    if not m:
        return None
    parts = m.group(1).strip().split()
    return [float(x) for x in parts]

def extract_coeffs(path, species):
    with open(path, 'r') as f:
        lines = f.readlines()
    i = 0
    n = len(lines)
    while i < n:
        if lines[i].strip() == species:
            # find '{'
            i += 1
            while i < n and '{' not in lines[i]:
                i += 1
            i += 1
            low = high = None
            while i < n and '}' not in lines[i]:
                ln = lines[i].strip()
                if ln.startswith("kappaLowT"):
                    low = parse_paren_nums(ln)
                elif ln.startswith("kappaHighT"):
                    high = parse_paren_nums(ln)
                i += 1
            return low, high
        i += 1
    return None, None

def eval_poly_asc(coeffs, T):
    # coeffs: [a0, a1, ...] ascending a0 + a1*T + a2*T^2 ...
    coeffs = np.asarray(coeffs, dtype=float)
    T = np.array(T, dtype=float)
    val = np.zeros_like(T)
    for i, a in enumerate(coeffs):
        val += a * (T**i)
    return val

def compute_cantera_kappa(gas, species, T):
    vals = np.empty_like(T)
    gas.X = species + ":1.0"
    gas.transport_model = "mixture-averaged"
    for i, TT in enumerate(T):
        gas.TP = TT, ct.one_atm
        vals[i] = gas.thermal_conductivity
    return vals

def plot_species_with_cantera(ax, species, low_coeffs, high_coeffs, gas):
    T_full = np.linspace(200, 4000, 1000)
    # Cantera reference (compute once)
    can_vals = compute_cantera_kappa(gas, species, T_full)

    # plot Cantera as points (small)
    ax.plot(T_full, can_vals, color='k', marker='.', linestyle='None', markersize=3, label=f"{species} Cantera")

    # plot fitted polynomials only on their intervals (solid)
    if low_coeffs:
        T_low = T_full[T_full <= SPLIT_T]
        y_low = eval_poly_asc(low_coeffs, T_low)
        y_low = np.where(np.isfinite(y_low), y_low, np.nan)
        ax.plot(T_low, y_low, '-', linewidth=1.5, label=f"{species} fit low")
    if high_coeffs:
        T_high = T_full[T_full >= SPLIT_T]
        y_high = eval_poly_asc(high_coeffs, T_high)
        y_high = np.where(np.isfinite(y_high), y_high, np.nan)
        ax.plot(T_high, y_high, '-', linewidth=1.5, label=f"{species} fit high")

    ax.axvline(SPLIT_T, color='gray', linestyle='--', linewidth=0.8)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Thermal conductivity (W/m/K)")
    ax.set_title(species)
    ax.legend(fontsize='small')
    ax.grid(True)
    # keep reasonable y-limits: center on Cantera data to avoid explosion by polynomial outside interval
    finite_can = can_vals[np.isfinite(can_vals)]
    if finite_can.size:
        ymin = np.min(finite_can)*0.5
        ymax = np.max(finite_can)*1.5
        if ymin < ymax:
            ax.set_ylim(ymin, ymax)

def main():
    try:
        gas = ct.Solution(MECH)
    except Exception as e:
        print("加载 Cantera 机制失败:", e)
        return

    fig, axes = plt.subplots(1, len(SPECIES), figsize=(6*len(SPECIES), 4))
    if len(SPECIES) == 1:
        axes = [axes]
    for ax, sp in zip(axes, SPECIES):
        low, high = extract_coeffs(THERMO_FILE, sp)
        if low is None and high is None:
            print(f"{sp}: 未找到系数 in {THERMO_FILE}")
            continue
        print(f"{sp}: kappaLowT len={len(low) if low else None}, kappaHighT len={len(high) if high else None}")
        plot_species_with_cantera(ax, sp, low, high, gas)

    plt.tight_layout()
    plt.savefig("kappa_O2_N2_with_cantera.png", dpi=200)
    print("图像已保存为 kappa_O2_N2_with_cantera.png")
    plt.show()

if __name__ == "__main__":
    main()