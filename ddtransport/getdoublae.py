import cantera as ct 
import numpy as np

# ===== Setting ====== #
mech = "chem.yaml"
Tlow = 200
Thigh = 6000
grids = 200

# 分段拟合开关及区间设置
split_kappa = True  # True=分段，False=单区间
low_T_range = (200, 1000)
high_T_range = (1000, 6000)
low_grids = 100
high_grids = 100

def transportTemplate(species_name, muCoeffs_arr, kappaLow_arr=None, kappaHigh_arr=None):
    def fmt_coeffs(arr):
        return " ".join(f"{v:.5e}" for v in arr)

    mu_block = "        muCoeffs<8>       (" + fmt_coeffs(muCoeffs_arr) + ");\n"

    if (kappaLow_arr is None) and (kappaHigh_arr is None):
        kappa_block = ""
    elif kappaHigh_arr is None:
        # single interval: write as kappaCoeffs (ascending order expected by reader)
        kappa_block = "        kappaCoeffs<8>    (" + fmt_coeffs(kappaLow_arr) + ");\n"
    else:
        # split: write low then high (both ascending a0..a7)
        kappa_block = (
            "        kappaLowT<8>    (" + fmt_coeffs(kappaLow_arr) + ");\n"
            "        kappaHighT<8>   (" + fmt_coeffs(kappaHigh_arr) + ");\n"
        )

    string = (
        f"{species_name}\n"
        "{\n"
        "    transport\n"
        "    {\n"
        f"{mu_block}"
        f"{kappa_block}"
        "    }\n"
        "}\n"
    )
    return string

with open("thermoPolynomial","w+") as f:
    gas = ct.Solution(mech)
    for species in gas.species_names:
        reactants_species = "{}:1".format(species)
        # 粘度全区间，用原始 coeff 数组（降幂）
        TList = np.linspace(Tlow, Thigh, grids)
        muList = []
        kappaList = []
        for T in TList:
            gas.TPY = T, 101325, reactants_species
            gas.transport_model = "mixture-averaged"
            mu = gas.viscosity
            kappa = gas.thermal_conductivity
            muList.append(mu)
            kappaList.append(kappa)
        pmu = np.polyfit(TList, muList, 7)            # pmu: [a7..a0]
        # For mu we will write ascending order, so reverse
        pmu_asc = pmu[::-1]                            # [a0..a7]

        pmu_poly = np.poly1d(pmu)
        mu_rmse = np.sqrt(np.mean((pmu_poly(TList) - muList)**2))
        mu_max_rel = np.max(np.abs((pmu_poly(TList) - muList)/muList))

        if not split_kappa:
            pkappa = np.polyfit(TList, kappaList, 7)  # pkappa: [a7..a0]
            pkappa_asc = pkappa[::-1]                 # [a0..a7] for file
            pkappa_poly = np.poly1d(pkappa)
            kappa_rmse = np.sqrt(np.mean((pkappa_poly(TList) - kappaList)**2))
            kappa_max_rel = np.max(np.abs((pkappa_poly(TList) - kappaList)/kappaList))
            print(f"{species}: mu_RMSE={mu_rmse:.3e}, mu_max_rel={mu_max_rel:.2%}; "
                  f"kappa_RMSE={kappa_rmse:.3e}, kappa_max_rel={kappa_max_rel:.2%}")
            f.write(transportTemplate(species, pmu_asc, pkappa_asc, None))
        else:
            # 低温区（返回降幂数组）
            TList_low = np.linspace(low_T_range[0], low_T_range[1], low_grids)
            kappaList_low = []
            for T in TList_low:
                gas.TPY = T, 101325, reactants_species
                gas.transport_model = "mixture-averaged"
                kappa = gas.thermal_conductivity
                kappaList_low.append(kappa)
            pkappa_low = np.polyfit(TList_low, kappaList_low, 7)   # [a7..a0]
            pkappa_low_poly = np.poly1d(pkappa_low)
            kappa_rmse_low = np.sqrt(np.mean((pkappa_low_poly(TList_low) - kappaList_low)**2))
            kappa_max_rel_low = np.max(np.abs((pkappa_low_poly(TList_low) - kappaList_low)/np.array(kappaList_low)))
            pkappa_low_asc = pkappa_low[::-1]                     # [a0..a7]

            # 高温区
            TList_high = np.linspace(high_T_range[0], high_T_range[1], high_grids)
            kappaList_high = []
            for T in TList_high:
                gas.TPY = T, 101325, reactants_species
                gas.transport_model = "mixture-averaged"
                kappa = gas.thermal_conductivity
                kappaList_high.append(kappa)
            pkappa_high = np.polyfit(TList_high, kappaList_high, 7)  # [a7..a0]

            # 保证在交界点连续（C0），注意 pkappa_low/pkappa_high 都是降幂数组
            Tj = low_T_range[1]
            val_low_at_Tj = np.polyval(pkappa_low, Tj)
            val_high_at_Tj = np.polyval(pkappa_high, Tj)
            diff = val_low_at_Tj - val_high_at_Tj
            pkappa_high[-1] += diff   # adjust constant term (降幂数组最后一项为 a0)
            pkappa_high_poly = np.poly1d(pkappa_high)
            kappa_rmse_high = np.sqrt(np.mean((pkappa_high_poly(TList_high) - kappaList_high)**2))
            denom_high = np.maximum(np.abs(kappaList_high), 1e-30)
            kappa_max_rel_high = np.max(np.abs((pkappa_high_poly(TList_high) - kappaList_high)/denom_high))

            # convert to ascending (a0..a7) for writing to file
            pkappa_high_asc = pkappa_high[::-1]

            print(f"{species}: mu_RMSE={mu_rmse:.3e}, mu_max_rel={mu_max_rel:.2%}; "
                  f"kappa_low_RMSE={kappa_rmse_low:.3e}, kappa_low_max_rel={kappa_max_rel_low:.2%}; "
                  f"kappa_high_RMSE={kappa_rmse_high:.3e}, kappa_high_max_rel={kappa_max_rel_high:.2%}")

            # write with correct order: mu (asc), kappaLow (asc), kappaHigh (asc)
            f.write(transportTemplate(species, pmu_asc, pkappa_low_asc, pkappa_high_asc))