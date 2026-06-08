"""
MAMC: Applegate Mechanism & Lanza Spin-Orbit Coupling Calculator
===============================================================

This core module implements and evaluates the energy requirements of various 
magnetic quadrupole moment and spin-orbit coupling variations observed in 
Post-Common Envelope Binary (PCEB) systems.

It follows the numerical frameworks established by:
  - Tian et al. (2009)
  - Völschow et al. (2016)
  - Lanza (2020)
"""

import math
import numpy as np
from astropy.constants import G, sigma_sb, M_sun, R_sun, L_sun
from astropy import units as u

# Global physical constants converted to CGS units
G_cgs = G.cgs.value
sigma_SB = sigma_sb.cgs.value
M_sun_cgs = M_sun.cgs.value
R_sun_cgs = R_sun.cgs.value
L_sun_cgs = L_sun.cgs.value

# Time unit conversions
day_to_sec = u.day.to(u.s)
year_to_sec = u.year.to(u.s)
year_to_day = u.year.to(u.d)

# IAU 2015 Resolution B3 nominal solar effective temperature
T_sun = 5772  


class ApplegateMechanism:
    """
    Encapsulates analytical formulations for evaluating the structural and energetic 
    feasibility of the Applegate mechanism in close binary systems.
    """

    @staticmethod
    def _resolve_dP_over_P(P_mod_yr: float, A_sec: float = None, dP_over_P: float = None) -> float:
        """
        Unifies amplitude tracking configurations into a single dimensionless 
        fractional period variation ratio.

        Args:
            P_mod_yr (float): The long-term period modulation cycle length in years.
            A_sec (float, optional): The timing semi-amplitude in seconds.
            dP_over_P (float, optional): Direct fractional period change ratio.

        Returns:
            float: The full-amplitude relative period variation ratio over the cycle.

        Raises:
            ValueError: If neither `A_sec` nor `dP_over_P` is supplied.
        """
        if A_sec is None and dP_over_P is None:
            raise ValueError("Amplitude tracking missing! Provide 'A_sec' OR 'dP_over_P'.")
        if dP_over_P is not None:
            return float(dP_over_P)
        
        # Scale timing semi-amplitude to full amplitude variation baseline
        return 2.0 * math.pi * (A_sec * 2) / (P_mod_yr * year_to_sec)

    @classmethod
    def thin_shell_ratio(cls, M_sec: float, R_sec: float, T_sec: float, a_bin: float, 
                         P_bin: float, P_mod_yr: float, A_sec: float = None, dP_over_P: float = None) -> float:
        """
        Computes the energy constraint ratio using the thin-shell approximation 
        derived by Tian et al. (2009).

        This model evaluates the absolute period change against the total energy 
        available from the secondary component over one full modulation cycle.
        
        .. math::
           \\frac{\\Delta E}{E_{\\mathrm{sec}}} = 0.233 \\left( \\frac{M_{\\mathrm{sec}}}{M_{\\odot}} \\right)^{3}
           \\left( \\frac{R_{\\mathrm{sec}}}{R_{\\odot}} \\right)^{-10} \\left( \\frac{T_{\\mathrm{sec}}}{6000~\\mathrm{K}} \\right)^{-4}
           \\left( \\frac{a_{\\mathrm{bin}}}{R_{\\odot}} \\right)^{4} \\left( \\frac{\\Delta P}{\\mathrm{s}} \\right)^{2}
           \\left( \\frac{P_{\\mathrm{mod}}}{\\mathrm{yr}} \\right)^{-1}

        Args:
            M_sec (float): Mass of the secondary star in solar units (:math:`M_\\odot`).
            R_sec (float): Radius of the secondary star in solar units (:math:`R_\\odot`).
            T_sec (float): Effective temperature of the active secondary star in Kelvin.
            a_bin (float): Binary separation distance in solar radii (:math:`R_\\odot`).
            P_bin (float): Binary orbital period in days.
            P_mod_yr (float): Long-term period modulation cycle length in years.
            A_sec (float, optional): Timing semi-amplitude in seconds.
            dP_over_P (float, optional): Direct relative period change.

        Returns:
            float: Dimensionless energy ratio (:math:`\\Delta E / E_{\\mathrm{sec}}`). Values >> 1 imply physical unfeasibility.
        """
        internal_dP = cls._resolve_dP_over_P(P_mod_yr, A_sec, dP_over_P)
        delta_P_sec = internal_dP * P_bin * day_to_sec
        
        return (0.233 * (M_sec**3) * (R_sec**-10) * (T_sec / 6000)**(-4) * (a_bin**4) * (delta_P_sec**2) * (P_mod_yr**-1))

    @classmethod
    def constant_density_ratio(cls, M_sec: float, R_sec: float, L_sec: float, a_bin: float, 
                               P_mod_yr: float, A_sec: float = None, dP_over_P: float = None) -> float:
        """
        Computes the energy validation threshold assuming a constant-density stellar profile.

        Applies the analytical scaling relation derived in Eq. (22) of Völschow et al. (2016).
        
        .. math::
           \\frac{\\Delta E}{E_{\\mathrm{sec}}} \\simeq 1.1 \\times 10^{7} \\left( \\frac{\\Delta P}{P_{\\mathrm{bin}}} \\right)
           \\left( \\frac{a_{\\mathrm{bin}}}{R_{\\odot}} \\right)^{2} \\left( \\frac{M_{\\mathrm{sec}}}{M_{\\odot}} \\right)^{2}
           \\left( \\frac{R_{\\mathrm{sec}}}{R_{\\odot}} \\right)^{-3} \\left( \\frac{P_{\\mathrm{mod}}}{\\mathrm{yr}} \\right)^{-1}
           \\left( \\frac{L_{\\mathrm{sec}}}{L_{\\odot}} \\right)^{-1}

        Args:
            M_sec (float): Mass of the secondary star in solar units (:math:`M_\\odot`).
            R_sec (float): Radius of the secondary star in solar units (:math:`R_\\odot`).
            L_sec (float): Luminosity of the active component in solar units (:math:`L_\\odot`).
            a_bin (float): Orbital separation distance in solar radii (:math:`R_\\odot`).
            P_mod_yr (float): Long-term period modulation cycle length in years.
            A_sec (float, optional): Timing semi-amplitude in seconds.
            dP_over_P (float, optional): Direct relative period change.

        Returns:
            float: Required relative energy efficiency budget ratio (:math:`\\Delta E / E_{\\mathrm{sec}}`).
        """
        internal_dP = cls._resolve_dP_over_P(P_mod_yr, A_sec, dP_over_P)
        
        return (1.1e7 * internal_dP * (a_bin**2) * (M_sec**2) * (R_sec**-3) * (P_mod_yr**-1) * (L_sec**-1))

    @classmethod
    def two_zone_ratio(cls, M_sec: float, R_sec: float, L_sec: float, a_bin: float, 
                       P_bin: float, P_mod_yr: float, A_sec: float = None, dP_over_P: float = None) -> tuple:
        """
        Evaluates system stability through the advanced two-zone variable shell core model.

        Models the active secondary as an inner dense core and an outer fluid convective envelope 
        exchanging angular momentum, changing the stellar quadrupole moment (Völschow et al. 2016).
        
        .. math::
           \\frac{\\Delta E}{E_{\\mathrm{sec}}} = k_{1}\\, \\frac{M_{\\mathrm{sec}} R_{\\mathrm{sec}}^{2}}{P_{\\mathrm{bin}}^{2}\\, P_{\\mathrm{mod}}\\, L_{\\mathrm{sec}}}
           \\left( 1 \\pm \\sqrt{1 - k_{2}\\, G\\, \\frac{a_{\\mathrm{bin}}^{2}\\, M_{\\mathrm{sec}}\\, P_{\\mathrm{bin}}^{2}}{R_{\\mathrm{sec}}^{5}} \\frac{\\Delta P}{P_{\\mathrm{bin}}}} \\right)^{2}

        Args:
            M_sec (float): Mass of the secondary star in solar units (:math:`M_\\odot`).
            R_sec (float): Radius of the secondary star in solar units (:math:`R_\\odot`).
            L_sec (float): Luminosity of the active component in solar units (:math:`L_\\odot`).
            a_bin (float): Semi-major axis separation distance in solar radii (:math:`R_\\odot`).
            P_bin (float): Orbital period of the close binary system in days.
            P_mod_yr (float): Long-term period modulation cycle length in years.
            A_sec (float, optional): Timing semi-amplitude in seconds.
            dP_over_P (float, optional): Direct relative period change.

        Returns:
            tuple:
                - **ratio** (*float*): Calculated energy efficiency ratio, or -1 if unfeasible.
                - **applegate_A** (*float*): Dimensionless structural coupling parameter :math:`A`. Physical solutions require :math:`A \\le 1`.
        """
        internal_dP = cls._resolve_dP_over_P(P_mod_yr, A_sec, dP_over_P)
        
        # Scaling inputs cleanly to absolute CGS units
        a_cgs = a_bin * R_sun_cgs
        M_cgs = M_sec * M_sun_cgs
        R_cgs = R_sec * R_sun_cgs
        L_cgs = L_sec * L_sun_cgs
        P_bin_sec = P_bin * day_to_sec
        P_mod_sec = P_mod_yr * year_to_day * day_to_sec

        # Empirically derived internal structure constants for low-mass MS stars
        k1, k2 = 0.133, 3.42
        
        applegate_A = k2 * G_cgs * (a_cgs**2 * P_bin_sec**2 * M_cgs / R_cgs**5) * internal_dP

        if applegate_A > 1.0:
            return -1, applegate_A

        prefactor = k1 * (M_cgs * R_cgs**2) / (P_bin_sec**2 * P_mod_sec * L_cgs)
        ratio = prefactor * (1.0 - math.sqrt(1.0 - applegate_A))**2

        return ratio, applegate_A


class LanzaMechanism:
    """
    Evaluates radiant energy constraints and mechanical shearing indicators 
    for the Lanza (2020) dynamic spin-orbit envelope coupling mechanism.
    """

    @classmethod
    def analyze(cls, P_bin: float, P_mod_yr: float, *, A_sec: float = None, dP_over_P: float = None, 
                M_pri: float = None, M_sec: float = None, R_sec: float = None, T_sec: float = None, 
                L_sec: float = None, a_bin: float = None, use_r0: str = 'a', k2: float = 0.2, 
                Ip_kgm2: float = None, mr0_sq_kgm2: float = None, mr0_Ip_ratio: float = None) -> tuple:
        """
        Executes a thermodynamic validation test for Lanza spin-orbit coupling.
        
        Maps observed tracking timing deviations down to the mechanical kinetic work done by 
        cyclic envelope speed modulations (:math:`\\Delta \\Omega / \\Omega`), verifying boundaries 
        against solar scaling parameters or tabular literature parameters.

        Args:
            P_bin (float): Binary system orbital period in days.
            P_mod_yr (float): Long-term period modulation sequence length in years.
            A_sec (float, optional): Observed timing semi-amplitude in seconds.
            dP_over_P (float, optional): Pre-computed relative period variation index.
            M_pri (float, optional): White dwarf primary component mass in solar units (:math:`M_\\odot`).
            M_sec (float, optional): Active low-mass secondary companion mass in solar units (:math:`M_\\odot`).
            R_sec (float, optional): Active companion radius in solar units (:math:`R_\\odot`).
            T_sec (float, optional): Active companion effective temperature in Kelvin.
            L_sec (float, optional): Active companion luminosity in solar units (:math:`L_\\odot`).
            a_bin (float, optional): Separational semi-major axis distance in solar radii (:math:`R_\\odot`).
            use_r0 (str, optional): Geometric constraint rule ('a' or 'a1'). Defaults to 'a'.
            k2 (float, optional): Gyration constant for the active convective shell. Defaults to 0.2.
            Ip_kgm2 (float, optional): Direct envelope moment of inertia invariant in :math:`\\text{kg}\\cdot\\text{m}^2`.
            mr0_sq_kgm2 (float, optional): Direct mass-radius framework scale in :math:`\\text{kg}\\cdot\\text{m}^2`.
            mr0_Ip_ratio (float, optional): Pre-computed dimensionless coupling ratio metric.

        Returns:
            tuple:
                - **deltaE_over_L_yr** (*float*): Energy ratio normalized across 1 Earth year.
                - **deltaE_over_LPmod** (*float*): Energy ratio normalized across full modulation cycle timeline.
                - **dOmega_over_Omega** (*float*): Internal required mechanical kinetic velocity shift.
                - **internal_dP** (*float*): Final active relative variation baseline.

        Raises:
            ValueError: If inputs break structural parameter dependency groupings.
        """
        # 1. Validation Guard Clauses for Mutually Exclusive Options
        if A_sec is None and dP_over_P is None:
            raise ValueError("Variation trace missing! Supply 'A_sec' OR 'dP_over_P'.")
            
        if Ip_kgm2 is None:
            if None in (M_pri, M_sec, R_sec, a_bin):
                raise ValueError("Observational parameter mismatch! Supply ('M_pri', 'M_sec', 'R_sec', 'a_bin') when 'Ip_kgm2' is absent.")
        else:
            if mr0_Ip_ratio is None and mr0_sq_kgm2 is None:
                raise ValueError("Literature shortcuts incomplete! Supply 'mr0_Ip_ratio' OR 'mr0_sq_kgm2' when 'Ip_kgm2' is used.")

        if L_sec is None and (R_sec is None or T_sec is None):
            raise ValueError("Luminosity metrics missing! Provide 'L_sec' directly OR pass 'R_sec' and 'T_sec' to scale.")

        if use_r0.lower() not in ("a", "a1"):
            raise ValueError("Parameter 'use_r0' constraint mode must be either 'a' or 'a1'.")

        # 2. Time & Amplitude Scalings
        P_sec = P_bin * day_to_sec
        Pmod_sec = P_mod_yr * year_to_day * day_to_sec
        
        internal_dP = float(dP_over_P) if dP_over_P is not None else 2.0 * np.pi * (A_sec * 2) / Pmod_sec

        # 3. Luminosity Extraction
        if L_sec is not None:
            L_sec_cgs = L_sec * L_sun_cgs
        else:
            L_sec_solar = (R_sec ** 2) * ((T_sec / T_sun) ** 4)
            L_sec_cgs = L_sec_solar * L_sun_cgs

        # 4. Moment of Inertia (Ip) Configuration
        if Ip_kgm2 is not None:
            Ip = float(Ip_kgm2) * 1e7  
        else:
            Ip = k2 * (M_sec * M_sun_cgs) * ((R_sec * R_sun_cgs) ** 2)

        # 5. Spin-Orbit Coupling Scale Factor Matching
        if mr0_Ip_ratio is not None:
            ratio_mr0_over_Ip = float(mr0_Ip_ratio)
        elif mr0_sq_kgm2 is not None:
            ratio_mr0_over_Ip = (float(mr0_sq_kgm2) * 1e7) / Ip
        else:
            M_sec_cgs = M_sec * M_sun_cgs
            M_pri_cgs = M_pri * M_sun_cgs
            a_bin_cgs = a_bin * R_sun_cgs
            
            mu = (M_sec_cgs * M_pri_cgs) / (M_sec_cgs + M_pri_cgs)
            r0 = a_bin_cgs if use_r0.lower() == "a" else a_bin_cgs * (M_pri_cgs / (M_sec_cgs + M_pri_cgs))
            ratio_mr0_over_Ip = (mu * (r0 ** 2)) / Ip

        # 6. Shear & Normalized Radiant Energy Feasibility Indexes
        dOmega_over_Omega = - (ratio_mr0_over_Ip / 3.0) * internal_dP
        Omega = 2.0 * np.pi / P_sec
        dE_erg = Ip * Omega * abs(Omega * dOmega_over_Omega)

        deltaE_over_L_yr = (dE_erg / L_sec_cgs) / year_to_sec
        deltaE_over_LPmod = deltaE_over_L_yr / P_mod_yr

        return deltaE_over_L_yr, deltaE_over_LPmod, dOmega_over_Omega, internal_dP

class AzimuthalDynamoWave:
    """
    Encapsulates the analytical scaling framework for the Azimuthal Dynamo Wave (ADW) 
    mechanism following Navarrete et al. (2026).
    
    This model evaluates whether observed non-axisymmetric quadrupole moment changes, 
    driven by strong azimuthally migrating magnetic fields in rapidly rotating, 
    fully convective stars ($\alpha^2$ dynamo), can reproduce observed Eclipse 
    Timing Variation (ETV) amplitudes.
    """

    @staticmethod
    def calculate_rotation_rate(R_sec: float, P_bin: float) -> tuple:
        """
        Computes the linear and relative angular rotation rates of the active secondary 
        star assuming strict tidal synchronization.

        Args:
            R_sec (float): Radius of the active secondary star in solar units ($R_\odot$).
            P_bin (float): Binary orbital period in days.

        Returns:
            tuple:
                - **omega_lin** (*float*): Linear equatorial rotational velocity in cm/s.
                - **omega_ratio** (*float*): Angular rotation velocity normalized to the solar rotation rate ($\Omega / \Omega_\odot$).
        """
        # Sidereal rotation period of the Sun in days (Snodgrass 1984)
        Prot_sun = 25.05 
        
        # Convert periods to absolute seconds
        P_bin_sec = P_bin * day_to_sec
        Prot2_sec = P_bin_sec  # Tidal locking assumption sets P_rot = P_bin
        
        # R_sec converted to centimeters
        R_sec_cm = R_sec * R_sun_cgs
        
        # Rotational velocities
        omega_lin = (2.0 * math.pi * R_sec_cm) / P_bin_sec
        omega_ratio = (2.0 * math.pi / Prot2_sec) / (2.0 * math.pi / (Prot_sun * day_to_sec))
        
        return omega_lin, omega_ratio

    @classmethod
    def analyze(cls, M_sec: float, R_sec: float, *, a_bin: float, P_bin: float, P_mod_yr: float, 
                f_scale: float = 1.0) -> tuple:
        """
        Predicts the theoretical O-C timing variation semi-amplitude generated by 
        non-axisymmetric mass deformations.
        
        Utilizes baseline principal inertia tensor variations ($\Delta I_{xx}, \Delta I_{yy}, \Delta I_{zz}$) 
        extracted from the fastest rotating 3D MHD model available ($10\Omega_\odot$) in 
        Navarrete et al. (2026).
        
        .. math::
           \Delta Q_{xx} = \Delta I_{xx} - \frac{1}{3}(\Delta I_{xx} + \Delta I_{yy} + \Delta I_{zz})
           
        .. math::
           O-C = \frac{9}{2\pi} \frac{\Delta Q_{xx}}{a^2 M_2} P_{\mathrm{mod}}

        Note:
            A negative $\Delta Q_{xx}$ is a real physical result of the dynamo simulations. 
            As rotation shifts to the rapid regime, strong magnetic fields concentrate 
            near the poles, inducing a sign flip that alters deformations from oblate 
            to prolate geometry relative to the rotational axis.

        Args:
            M_sec (float): Mass of the secondary star in solar units ($M_\odot$).
            R_sec (float): Radius of the secondary star in solar units ($R_\odot$).
            a_bin (float): Binary separation distance in solar radii ($R_\odot$).
            P_bin (float): Binary orbital period in days.
            P_mod_yr (float): Long-term period modulation cycle length in years.
            f_scale (float, optional): Dynamo configuration scaling coefficient parameter. 
                Defaults to 1.0 (absolute weakest configuration baseline).

        Returns:
            tuple:
                - **delta_Q_xx** (*float*): Non-axisymmetric quadrupole moment variation in absolute CGS ($g \cdot cm^2$).
                - **amp_OC** (*float*): Predicted theoretical O-C timing variation semi-amplitude in seconds.
                - **omega_lin_kms** (*float*): Linear tracking velocity of the active component in km/s.
                - **omega_ratio** (*float*): Relative angular rotation factor ($\Omega / \Omega_\odot$).
        """
        # Resolve tracking rotation dynamics
        omega_lin, omega_ratio = cls.calculate_rotation_rate(R_sec, P_bin)
        omega_lin_kms = omega_lin / 1e5  # cm/s to km/s
        
        # Absolute timeline conversions
        P_mod_sec = P_mod_yr * year_to_day * day_to_sec
        
        # Absolute metric conversions to CGS baselines
        a_cgs = a_bin * R_sun_cgs
        M_cgs = M_sec * M_sun_cgs

        # Baseline inertia tensor variations from Table 3 of Navarrete et al. (2026)
        # Tabulated in SI units (kg * m^2), converted to CGS (g * cm^2) by multiplying by 1e7
        delta_I_xx = 6.00192126e40 * f_scale * 1e7
        delta_I_yy = 6.00267489e40 * f_scale * 1e7
        delta_I_zz = 1.40925322e41 * f_scale * 1e7
        
        # Map structural deformations to non-axisymmetric quadrupole moment components
        delta_Q_xx = delta_I_xx - (delta_I_xx + delta_I_yy + delta_I_zz) / 3.0
        
        # Analytical approximation mapping mass shifts to timing amplitudes
        amp_OC = (9.0 / (2.0 * math.pi)) * ((delta_Q_xx * P_mod_sec) / ((a_cgs ** 2) * M_cgs))
        
        return delta_Q_xx, amp_OC, omega_lin_kms, omega_ratio
