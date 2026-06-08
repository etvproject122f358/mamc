"""
example_mamc: Comprehensive Diagnostics Pipeline for PCEB Timing Models
========================================================================

Main execution script for the MAMC framework. Processes data registries
and isolates outputs into structured tables with interpretation instructions.
"""

import pandas as pd
from target_params import TARGETS
from mamc import ApplegateMechanism, LanzaMechanism, AzimuthalDynamoWave

def run_diagnostics_pipeline():
    """Executes separated structural testing across the target catalog."""
    
    # Storage lists for isolated model metrics
    applegate_data = []
    lanza_data = []
    adw_data = []

    for name, sys in TARGETS.items():
        # ------------------------------------------------------------
        # 1. PROCESS APPLEGATE MECHANISM MODELS
        # ------------------------------------------------------------
        r_thin = ApplegateMechanism.thin_shell_ratio(
            M_sec=sys.M_sec, R_sec=sys.R_sec, T_sec=sys.T_sec, a_bin=sys.a_bin, 
            P_bin=sys.P_bin_days, P_mod_yr=sys.P_mod_yr, A_sec=sys.A_sec
        )
        r_const = ApplegateMechanism.constant_density_ratio(
            M_sec=sys.M_sec, R_sec=sys.R_sec, L_sec=sys.L_sec, a_bin=sys.a_bin, 
            P_mod_yr=sys.P_mod_yr, A_sec=sys.A_sec
        )
        r_two, app_A = ApplegateMechanism.two_zone_ratio(
            M_sec=sys.M_sec, R_sec=sys.R_sec, L_sec=sys.L_sec, a_bin=sys.a_bin, 
            P_bin=sys.P_bin_days, dP_over_P=sys.dP_over_Pbin, P_mod_yr=sys.P_mod_yr
        )
        
        applegate_data.append({
            "System": sys.name,
            "Applegate A": app_A,
            "Thin-Shell": r_thin,
            "Const-Density": r_const,
            "Two-Zone": r_two
        })

        # ------------------------------------------------------------
        # 2. PROCESS LANZA SPIN-ORBIT COUPLING MODEL
        # ------------------------------------------------------------
        l_yr, l_mod, dOm, _ = LanzaMechanism.analyze(
            P_bin=sys.P_bin_days, P_mod_yr=sys.P_mod_yr, A_sec=sys.A_sec,
            M_pri=sys.M_pri, M_sec=sys.M_sec, R_sec=sys.R_sec, T_sec=sys.T_sec, 
            L_sec=sys.L_sec, a_bin=sys.a_bin
        )
        
        lanza_data.append({
            "System": sys.name,
            "Observed A (s)": sys.A_sec,
            "dOmega/Omega": dOm,
            "dE/L (1 yr)": l_yr,
            "dE/L (Pmod)": l_mod
        })

        # ------------------------------------------------------------
        # 3. PROCESS AZIMUTHAL DYNAMO WAVE MODEL
        # ------------------------------------------------------------
        dQ_xx, theoretical_OC, v_lin, v_ratio = AzimuthalDynamoWave.analyze(
            M_sec=sys.M_sec, R_sec=sys.R_sec, a_bin=sys.a_bin, 
            P_bin=sys.P_bin_days, P_mod_yr=sys.P_mod_yr, f_scale=1.0
        )
        
        adw_data.append({
            "System": sys.name,
            "Observed A (s)": sys.A_sec,
            "V_eq (km/s)": v_lin,
            "Omega/Omega_sun": v_ratio,
            "Delta Q_xx": dQ_xx,
            "Expected O-C (s)": theoretical_OC
        })

    # Global display settings for standard scannability
    pd.options.display.float_format = "{:.3e}".format
    
    # ----------------------------------------------------------------
    # RENDER SECTION 1: APPLEGATE MECHANISM BLOCK
    # ----------------------------------------------------------------
    print("\n" + "="*85)
    print("1. APPLEGATE MECHANISM ENERGY FEASIBILITY CRITERIA (Völschow et al. 2016)")
    print("="*85)
    df_app = pd.DataFrame(applegate_data)
    print(df_app.to_string(index=False, formatters={
        "Applegate A": "{:.3f}".format, "Thin-Shell": "{:.2f}".format,
        "Const-Density": "{:.2f}".format, "Two-Zone": "{:.2f}".format
    }))
    print("\n>>> INSTRUCTIONS FOR INTERPRETATION:")
    print("  • Applegate Parameter (A): Must be <= 1.0. If A > 1, the model is physically")
    print("    disproven (Two-Zone returns -1) as the orbit cannot support the deformation.")
    print("  • Energy Ratios (dE/E_sec): Measures required shell work vs. nuclear energy budget.")
    print("    - dE/E_sec << 1 : Highly viable framework (e.g., QS Vir).")
    print("    - dE/E_sec ~ 1  : Borderline case; energetically tight but plausible.")
    print("    - dE/E_sec >> 1 : Disproven on energetic grounds; requires more energy than the")
    print("                      star produces over the cycle (e.g., NN Ser, HW Vir).")

    # ----------------------------------------------------------------
    # RENDER SECTION 2: LANZA SPIN-ORBIT COUPLING BLOCK
    # ----------------------------------------------------------------
    print("\n" + "="*85)
    print("2. LANZA SPIN-ORBIT COUPLING THERMODYNAMIC BUDGETS (Lanza 2020)")
    print("="*85)
    df_lan = pd.DataFrame(lanza_data)
    print(df_lan.to_string(index=False, formatters={
        "Observed A (s)": "{:.2f}".format, "dOmega/Omega": "{:.3e}".format,
        "dE/L (1 yr)": "{:.2f}".format, "dE/L (Pmod)": "{:.2f}".format
    }))
    print("\n>>> INSTRUCTIONS FOR INTERPRETATION:")
    print("  • dE/L (Pmod): The fraction of total radiant energy that must be channeled into")
    print("    mechanical envelope shear over the modulation period.")
    print("    - dE/L (Pmod) <= 0.001 : Highly feasible; requires < 0.1% of radiant energy.")
    print("    - dE/L (Pmod) ~ 0.01-0.1: Demanding; requires 1% to 10% of total stellar luminosity.")
    print("    - dE/L (Pmod) >= 1.0   : Physically impossible; requires converting more energy than")
    print("                             the active core generates into pure kinetic work.")

    # ----------------------------------------------------------------
    # RENDER SECTION 3: AZIMUTHAL DYNAMO WAVE BLOCK
    # ----------------------------------------------------------------
    print("\n" + "="*85)
    print("3. NON-AXISYMMETRIC AZIMUTHAL DYNAMO WAVE AMPLITUDES (Navarrete et al. 2026)")
    print("="*85)
    df_adw = pd.DataFrame(adw_data)
    print(df_adw.to_string(index=False, formatters={
        "Observed A (s)": "{:.2f}".format, "V_eq (km/s)": "{:.2f}".format,
        "Omega/Omega_sun": "{:.2f}".format, "Delta Q_xx": "{:.2e}".format,
        "Expected O-C (s)": "{:.2f}".format
    }))
    print("\n>>> INSTRUCTIONS FOR INTERPRETATION:")
    print("  • Delta Q_xx: Calculated non-axisymmetric quadrupole moment change in the orbital plane.")
    print("    - Negative sign is a real physical trait indicating that rapid rotators switch from")
    print("      oblate (pancake) to prolate (cigar-shaped) polar magnetic concentrations.")
    print("  • Expected O-C (s) vs Observed A (s): Evaluation of wave amplitude driving potential.")
    print("    - This run fixes f_scale = 1.0, modeling the absolute weakest baseline dynamo.")
    print("    - If Expected O-C meets or exceeds the Observed A baseline even under f=1.0, the ADW")
    print("      mechanism stands as a robust, viable driver for the observed ETV curves.")
    print("="*85 + "\n")

if __name__ == "__main__":
    run_diagnostics_pipeline()
