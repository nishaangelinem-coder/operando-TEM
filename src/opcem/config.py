"""Central configuration for the OPCEM digital twin and analysis pipeline.

All physical parameters live here, not in analysis scripts (design rule 5 of
docs/04_experimental_setup.md).  Values are documented in
docs/04_experimental_setup.md section 5.4, including which are anchored to
published measurements and which are plausible simulation choices.

IMPORTANT: every number produced by this package is a SIMULATION output.
Nothing here is an experimental measurement.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Dict, Tuple

KB_EV = 8.617333262e-5  # eV/K
MASTER_SEED = 20260922


@dataclass
class SupportConfig:
    """Support lattice and site-type composition."""
    name: str = "CeO2(100)"
    lattice_nm: float = 0.383           # surface site spacing
    field_nm: float = 40.0              # imaged field of view edge length
    slab_thickness_nm: float = 8.0      # sets the depth range atoms occupy
    frac_vacancy: float = 0.06          # O-vacancy site fraction
    frac_step: float = 0.03             # step-edge site fraction
    # Extra diffusion barrier contributed by each site type (eV)
    e_trap: Dict[str, float] = field(default_factory=lambda: {
        "terrace": 0.00, "vacancy": 0.15, "step": 0.08})
    # Surface corrugation: z = amp * sin(2 pi x / wav) * sin(2 pi y / wav).
    # Lateral diffusion therefore carries a real out-of-plane component, which
    # projection along the optic axis destroys (this is what makes G3 testable).
    corrug_amp_nm: float = 0.70
    corrug_wavelength_nm: float = 5.0


@dataclass
class KineticsConfig:
    """Ground-truth Pt mobility / aggregation kinetics."""
    nu0_hop: float = 1.0e12             # s^-1 attempt frequency
    e_hop_terrace: float = 1.32         # eV, monomer hop barrier on terrace
    capture_radius_nm: float = 0.55     # encounter radius for Pt1 + Pt1
    # Binding energy holding the outermost atom in a Pt_n entity (eV)
    e_bind: Dict[int, float] = field(default_factory=lambda: {
        2: 1.46, 3: 1.60, 4: 1.70, 5: 1.76})
    e_bind_large: float = 1.82
    nu0_diss: float = 1.0e12
    # Oxidative redispersion: k_diss multiplied by (1 + gamma * sqrt(p_O2/bar))
    gamma_ox: float = 6.0
    # Beam contribution to hop rate: k_beam = beta_beam * dose_rate
    # dose_rate in e-/A^2/s ; beta chosen so that phi=1e4 gives a hop rate
    # comparable to the 300 C thermal rate (a deliberately strong confound).
    beta_beam: float = 2.0e-4


@dataclass
class ActivityConfig:
    """Ground-truth intrinsic activity of each Pt_n motif.

    The dimer is the most productive motif per Pt atom.  The pipeline is never
    told this; recovering it is the validation target.
    """
    # Pre-exponential per entity (molecules/s), indexed by n
    a_n: Dict[int, float] = field(default_factory=lambda: {
        1: 4.0e6, 2: 9.0e4, 3: 1.7e5, 4: 3.6e5})
    a_n_large_per_atom: float = 1.17e5  # per Pt atom for n >= 5
    e_act: Dict[int, float] = field(default_factory=lambda: {
        1: 0.82, 2: 0.55, 3: 0.68, 4: 0.74})
    e_act_large: float = 0.78
    # Site-type multiplier on activity (vacancy-anchored Pt is more active)
    site_mult: Dict[str, float] = field(default_factory=lambda: {
        "terrace": 1.0, "vacancy": 1.45, "step": 1.15})
    # Langmuir adsorption (bar^-1 at 300 C, van't Hoff corrected internally)
    k_ads_co: float = 55.0
    k_ads_o2: float = 6.5
    dh_ads_co: float = -1.30            # eV
    dh_ads_o2: float = -0.95            # eV
    t_ref_k: float = 573.15


@dataclass
class ReactorConfig:
    """Nanoreactor geometry, operating point, and transport field."""
    p_total_mbar: float = 700.0
    x_co: float = 0.02
    x_o2: float = 0.04
    t_setpoint_c: float = 300.0
    flow_sccm: float = 2.0
    # Imposed inlet->outlet gradients across the reactive zone
    t_gradient_c: float = 14.0          # outlet hotter than inlet
    product_enrichment: float = 0.18    # outlet enriched vs inlet (ref R4)
    zone_length_um: float = 200.0
    zone_width_um: float = 50.0
    # Ratio of the catalyst film's Pt-accessible surface area to the projected
    # area of the reactive zone.  A ~100 nm film of 5 nm CeO2 nanocrystals
    # presents of order 10^4-10^5 times the projected area (equivalently, a
    # ~7 ug, ~50 m^2/g powder over a 200 x 50 um zone).  The imaged field, by
    # contrast, is a single projected patch, so this factor enters the
    # representativeness fraction directly.
    catalyst_area_multiplier: float = 3.5e4
    # Representativeness: the imaged field's Pt-accessible area as a fraction
    # of the whole catalyst's Pt-accessible area.
    @property
    def projected_area_fraction(self) -> float:
        a_img = (40.0e-3) ** 2                        # um^2
        a_zone = self.zone_length_um * self.zone_width_um
        return a_img / a_zone

    @property
    def f_rep(self) -> float:
        return self.projected_area_fraction / self.catalyst_area_multiplier

    def molar_flow_mol_per_s(self) -> float:
        """Total molar flow. 1 sccm = 7.4362e-7 mol/s at 0 C, 1 atm."""
        return self.flow_sccm * 7.4362e-7


@dataclass
class TransportConfig:
    """Reactor+line+inlet impulse response h(t): gamma kernel + dead time."""
    tau0_s: float = 1.85                # pure dead time
    gamma_k: float = 2.6                # gamma shape
    gamma_theta_s: float = 0.52         # gamma scale
    # Flow scaling: tau scales as 1/flow
    ref_flow_sccm: float = 2.0

    def tau_dead(self) -> float:
        return self.tau0_s + self.gamma_k * self.gamma_theta_s

    def sigma_disp(self) -> float:
        return self.gamma_theta_s * self.gamma_k ** 0.5

    def scaled(self, flow_sccm: float) -> "TransportConfig":
        s = self.ref_flow_sccm / max(flow_sccm, 1e-6)
        return TransportConfig(tau0_s=self.tau0_s * s,
                               gamma_k=self.gamma_k,
                               gamma_theta_s=self.gamma_theta_s * s,
                               ref_flow_sccm=self.ref_flow_sccm)


@dataclass
class ImagingConfig:
    """Virtual ADF-STEM / EELS parameters."""
    dose_rate: float = 1.0e3            # e-/A^2/s
    fps: float = 5.0
    pixel_nm: float = 0.02          # 0.2 A sampling, required for single atoms
    render_field_nm: float = 8.0    # pixel-level rendering uses a crop
    probe_sigma_nm: float = 0.045
    ref_dose_rate: float = 1.0e3
    loc_sigma_nm_at_ref: float = 0.025  # 0.25 A localisation precision
    drift_nm_per_s: float = 0.004
    scan_jitter_nm: float = 0.010
    window_bg_level: float = 0.12
    # ADF contrast model: effective Pt cross-section into the ADF detector and
    # the fraction of the incident dose the SiN window + support scatter into
    # it.  Chosen so that a single Pt atom has CNR ~ 5 at the reference dose.
    sigma_adf_pt_A2: float = 0.125
    bg_scatter_frac: float = 0.05
    # Both calibrated against the pixel-level detector of vision.py over a
    # dose series (study S0): the CNR threshold reproduces the measured
    # monomer recall, and the false-positive density was measured on
    # atom-free frames (2.6e-3 nm^-2 at the reference dose).
    detect_threshold_cnr: float = 3.54
    false_positive_per_nm2: float = 2.6e-3
    fp_dose_exponent: float = 1.5
    depth_resolution_nm: float = 6.6    # single-projection ptychography (R8)
    snapshot3d_interval_s: float = 60.0
    eels_snr_at_ref: float = 18.0


@dataclass
class MSConfig:
    """Virtual quadrupole MS."""
    rate_hz: float = 10.0
    # Calibration factors: ion current per unit molar flow (arbitrary units)
    k_cal: Dict[str, float] = field(default_factory=lambda: {
        "CO": 1.00, "O2": 0.86, "CO2": 1.34, "Ar": 0.72})
    noise_rel: float = 0.020            # relative Gaussian noise
    drift_rel: float = 0.015            # 1/f baseline drift amplitude
    n2_interference: float = 0.004      # additive on m/z 28 channel
    # A quadrupole detects a CONCENTRATION, so its flux detection limit scales
    # with the total flow: flux_min = x_min * F.  1 ppm is optimistic for a
    # differentially pumped inlet with an electron multiplier.
    detection_mole_fraction: float = 1.0e-6
    background: Dict[str, float] = field(default_factory=lambda: {
        "CO": 0.010, "O2": 0.008, "CO2": 0.006, "Ar": 0.004})


@dataclass
class AnalysisConfig:
    """Pre-registered analysis thresholds."""
    bin_s: float = 1.0
    pre_window_s: float = 6.0
    post_window_s: float = 6.0
    event_guard_s: float = 2.0          # min separation between kept events
    n_bootstrap: int = 2000
    n_surrogate: int = 1000
    alpha: float = 0.05
    # correlate.py refuses to call a lag significant below this multiple of sigma_disp
    lag_significance_sigma_mult: float = 2.0
    bpi_exclusion_threshold: float = 0.20
    elastic_net_l1_ratio: float = 0.5
    hmm_n_states: int = 4
    tikhonov_lambda_grid: Tuple[float, ...] = (1e-4, 3e-4, 1e-3, 3e-3,
                                               1e-2, 3e-2, 1e-1, 3e-1, 1.0)


@dataclass
class ControlConfig:
    """Layer-3 policy parameters."""
    o2_pulse_duration_s: float = 4.0
    o2_pulse_x_o2: float = 0.20
    actuation_budget: int = 24          # max pulses per run (shared by all policies)
    lam_act: float = 0.0                # actuation cost weight (reported separately)
    trigger_n3_threshold: float = 1.5   # mean N>=3 entities over trigger window
    trigger_window_s: float = 10.0
    refractory_s: float = 20.0


@dataclass
class Config:
    support: SupportConfig = field(default_factory=SupportConfig)
    kinetics: KineticsConfig = field(default_factory=KineticsConfig)
    activity: ActivityConfig = field(default_factory=ActivityConfig)
    reactor: ReactorConfig = field(default_factory=ReactorConfig)
    transport: TransportConfig = field(default_factory=TransportConfig)
    imaging: ImagingConfig = field(default_factory=ImagingConfig)
    ms: MSConfig = field(default_factory=MSConfig)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)
    control: ControlConfig = field(default_factory=ControlConfig)
    n_pt_atoms: int = 90                # Pt atoms in the imaged field
    duration_s: float = 900.0
    seed: int = MASTER_SEED

    def as_dict(self) -> dict:
        return asdict(self)


DEFAULT = Config()

# Provenance tags enforced by descriptors.py (design rule 1)
PROVENANCE = {
    "N1": "census2d", "N2": "census2d", "N3": "census2d", "N4p": "census2d",
    "mean_size": "census2d", "cn_pt_pt": "census2d", "frac_dispersed": "census2d",
    "k_hop": "dynamic2d", "k_nuc": "dynamic2d", "k_diss": "dynamic2d",
    "tau_dimer": "dynamic2d", "msd_1s": "dynamic2d",
    "q_pt": "spectro", "theta_ads": "spectro",
    "T_c": "reactor", "p_co": "reactor", "p_o2": "reactor",
    "persistent_frac": "derived",
    "dose_rate": "reactor", "f_rep": "reactor",
    "cn_pt_support": "snapshot3d", "msd_3d_1s": "snapshot3d",
}
