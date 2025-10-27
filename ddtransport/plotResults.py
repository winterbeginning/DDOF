import re
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

THERMO_FILE = "thermoPolynomial"
MECH = "chem.yaml"
SPLIT_T = 1000.0

# 只画这些组分（可自定义）
SPECIES_LIST = [ "AL2O2", "AL2O3"]

def parse_paren_nums(s):
    m = re.search(r'\(([^)]*)\)', s)
    if not m:
        return None
    parts = m.group(1).strip().split()
    return [float(x) for x in parts]

def extract_all_species(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    species_blocks = {}
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i].strip()
        if line and not line.startswith("//") and not line.startswith("#") and not line.startswith("}") and not line.startswith("{"):
            species = line
            # find '{'
            i += 1
            while i < n and '{' not in lines[i]:
                i += 1
            i += 1
            block = []
            while i < n and '}' not in lines[i]:
                block.append(lines[i].strip())
                i += 1
            species_blocks[species] = block
        i += 1
    return species_blocks

def extract_coeffs(block):
    mu = kappa = kappa_low = kappa_high = None
    for ln in block:
        if ln.startswith("muCoeffs"):
            mu = parse_paren_nums(ln)
        elif ln.startswith("kappaCoeffs"):
            kappa = parse_paren_nums(ln)
        elif ln.startswith("kappaLowT"):
            kappa_low = parse_paren_nums(ln)
        elif ln.startswith("kappaHighT"):
            kappa_high = parse_paren_nums(ln)
    return mu, kappa, kappa_low, kappa_high

def eval_poly_asc(coeffs, T):
    coeffs = np.asarray(coeffs, dtype=float)
    T = np.array(T, dtype=float)
    val = np.zeros_like(T)
    for i, a in enumerate(coeffs):
        val += a * (T**i)
    return val

def compute_cantera(gas, species, T, prop):
    vals = np.empty_like(T)
    gas.X = species + ":1.0"
    gas.transport_model = "mixture-averaged"
    for i, TT in enumerate(T):
        gas.TP = TT, ct.one_atm
        if prop == "mu":
            vals[i] = gas.viscosity
        elif prop == "kappa":
            vals[i] = gas.thermal_conductivity
    return vals

def plot_species(ax, species, mu, kappa, kappa_low, kappa_high, gas):
    T_full = np.linspace(200, 4000, 1000)
    T_cantera = np.linspace(200, 4000, 50)

    # 粘度
    if mu:
        can_mu = compute_cantera(gas, species, T_cantera, "mu")
        y_mu = eval_poly_asc(mu, T_full)
        y_mu = np.where(np.isfinite(y_mu), y_mu, np.nan)
        ax[0].plot(T_full, y_mu, '-', color='b', linewidth=1.5, label=f"fit mu")
        ax[0].plot(T_cantera, can_mu, marker='o', ls='None', ms=4, mfc='none', mec='r', label="Cantera mu")
        ax[0].set_xlabel("Temperature (K)")
        ax[0].set_ylabel("Viscosity (Pa·s)")
        ax[0].set_title(f"{species} - Viscosity")
        ax[0].legend(fontsize='small')
        ax[0].grid(True)

    # 热导率
    if (kappa_low and kappa_high):
        can_kappa = compute_cantera(gas, species, T_cantera, "kappa")
        T_low = T_full[T_full <= SPLIT_T]
        T_high = T_full[T_full >= SPLIT_T]
        y_low = eval_poly_asc(kappa_low, T_low)
        y_high = eval_poly_asc(kappa_high, T_high)
        y_low = np.where(np.isfinite(y_low), y_low, np.nan)
        y_high = np.where(np.isfinite(y_high), y_high, np.nan)
        ax[1].plot(T_low, y_low, '-', color='b', linewidth=1.5, label="fit kappaLowT")
        ax[1].plot(T_high, y_high, '-', color='g', linewidth=1.5, label="fit kappaHighT")
        ax[1].axvline(SPLIT_T, color='gray', linestyle='--', linewidth=0.8)
        ax[1].plot(T_cantera, can_kappa, marker='o', ls='None', ms=4, mfc='none', mec='r', label="Cantera kappa")
        ax[1].set_xlabel("Temperature (K)")
        ax[1].set_ylabel("Thermal conductivity (W/m/K)")
        ax[1].set_title(f"{species} - Thermal conductivity")
        ax[1].legend(fontsize='small')
        ax[1].grid(True)
    elif kappa or kappa_low or kappa_high:
        # 单区间
        coeffs = kappa if kappa else (kappa_low if kappa_low else kappa_high)
        can_kappa = compute_cantera(gas, species, T_cantera, "kappa")
        y = eval_poly_asc(coeffs, T_full)
        y = np.where(np.isfinite(y), y, np.nan)
        ax[1].plot(T_full, y, '-', color='b', linewidth=1.5, label="fit kappa")
        ax[1].plot(T_cantera, can_kappa, marker='o', ls='None', ms=4, mfc='none', mec='r', label="Cantera kappa")
        ax[1].set_xlabel("Temperature (K)")
        ax[1].set_ylabel("Thermal conductivity (W/m/K)")
        ax[1].set_title(f"{species} - Thermal conductivity")
        ax[1].legend(fontsize='small')
        ax[1].grid(True)

def main():
    try:
        gas = ct.Solution(MECH)
    except Exception as e:
        print("加载 Cantera 机制失败:", e)
        return

    species_blocks = extract_all_species(THERMO_FILE)
    for sp in SPECIES_LIST:
        if sp not in species_blocks:
            print(f"{sp}: 未找到物性数据，跳过")
            continue
        block = species_blocks[sp]
        mu, kappa, kappa_low, kappa_high = extract_coeffs(block)
        print(f"{sp}: mu len={len(mu) if mu else None}, kappa len={len(kappa) if kappa else None}, "
              f"kappaLowT len={len(kappa_low) if kappa_low else None}, kappaHighT len={len(kappa_high) if kappa_high else None}")
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        plot_species(axes, sp, mu, kappa, kappa_low, kappa_high, gas)
        plt.tight_layout()
        plt.suptitle(sp)
        plt.show()

if __name__ == "__main__":
    main()