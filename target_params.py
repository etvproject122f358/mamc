"""
target_params: Post-Common Envelope Binary (PCEB) Data Registry
==============================================================

This registry holds historical, literature, and observational parameter sets 
for specific target eclipsing binary systems under tracking.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class BinarySystem:
    """
    Data structure defining empirical and orbital target attributes.
    
    Attributes:
        name (str): The common designated tracking name of the eclipsing binary system.
        M_pri (float): Mass of the white dwarf primary component in solar masses.
        M_sec (float): Mass of the active main-sequence companion in solar masses.
        R_sec (float): Radius of the active main-sequence companion in solar radii.
        T_sec (float): Effective blackbody temperature of the active companion in Kelvin.
        L_sec (float): Evaluated visual radiant luminosity profile in solar luminosities.
        a_bin (float): Calculated system separational distance scale in solar radii.
        P_bin_days (float): Empirical synodic orbital baseline frequency period in days.
        P_mod_yr (float): Cyclical timing variation anomaly wavelength period trace in years.
        dP_over_Pbin (float): Dimensionless fractional relative orbital period swing variant.
        A_sec (float): Extracted O-C observation diagram semi-amplitude depth in seconds.
        Ip_kgm2 (float, optional): Invariant shell inertia momentum scale listed in literature.
        mr0_Ip_ratio (float, optional): Invariant framework ratio metric listed in literature.
    """
    name: str
    M_pri: float
    M_sec: float
    R_sec: float
    T_sec: float
    L_sec: float
    a_bin: float
    P_bin_days: float
    P_mod_yr: float
    dP_over_Pbin: float
    A_sec: float
    Ip_kgm2: Optional[float] = None
    mr0_Ip_ratio: Optional[float] = None


# Populating catalog targets with verified observational records
TARGETS = {
    "NN Ser": BinarySystem(
        name="NN Ser", M_pri=0.535, M_sec=0.111, R_sec=0.149, T_sec=2920.0, L_sec=0.00147,
        a_bin=0.934, P_bin_days=0.130, P_mod_yr=15.482, dP_over_Pbin=7.1e-7, A_sec=27.65
    ),
    
    "HW Vir": BinarySystem(
        name="HW Vir", M_pri=0.485, M_sec=0.142, R_sec=0.175, T_sec=3084.0, L_sec=0.003,
        a_bin=0.860, P_bin_days=0.117, P_mod_yr=55.0, dP_over_Pbin=4.1e-6, A_sec=563.0
    ),
    
    "QS Vir": BinarySystem(
        name="QS Vir", M_pri=0.78, M_sec=0.43, R_sec=0.42, T_sec=3100.0, L_sec=0.0146,
        a_bin=1.27, P_bin_days=0.151, P_mod_yr=16.99, dP_over_Pbin=1.0e-6, A_sec=43.0
    ),
    
    "DD CrB": BinarySystem(
        name="DD CrB", M_pri=0.417, M_sec=0.127, R_sec=0.1619, T_sec=2357.0, L_sec=0.0007,
        a_bin=1.020, P_bin_days=0.16177, P_mod_yr=13.28487, dP_over_Pbin=7.055e-6, A_sec=5.10
    ),
    
    "NY Vir": BinarySystem(
        name="NY Vir", M_pri=0.471, M_sec=0.13, R_sec=0.155, T_sec=2048.0, L_sec=0.00038,
        a_bin=0.77, P_bin_days=0.1010159690, P_mod_yr=22.35, dP_over_Pbin=7.055e-6, A_sec=17.55
    )
}
