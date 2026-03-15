#!/usr/bin/env python3
"""
Thermodynamic analysis for CRISPR guide RNA validation.

This module implements thermodynamic calculations for gRNA-target binding,
including melting temperature, free energy, and specificity predictions.
Based on nearest-neighbor thermodynamic models for RNA-DNA hybrids.
"""

import numpy as np
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class ThermodynamicProfile:
    """Thermodynamic profile for a guide RNA-target interaction."""

    sequence: str
    target: str
    delta_g: float  # Free energy change (kcal/mol)
    delta_h: float  # Enthalpy change (kcal/mol)
    delta_s: float  # Entropy change (cal/mol·K)
    tm: float  # Melting temperature (°C)
    binding_efficiency: float  # Predicted binding efficiency (0-1)
    mismatch_penalty: float  # Free energy penalty from mismatches

    @property
    def stability_score(self) -> float:
        """Calculate overall stability score (0-1)."""
        # More negative ΔG = more stable
        # Normalize to 0-1 scale (-40 to 0 kcal/mol range)
        normalized_dg = max(0, min(1, (-self.delta_g) / 40))
        # High Tm is good (normalize 40-80°C to 0-1)
        normalized_tm = max(0, min(1, (self.tm - 40) / 40))
        # Combine with weights
        return 0.6 * normalized_dg + 0.4 * normalized_tm


class ThermodynamicCalculator:
    """Calculate thermodynamic parameters for oligonucleotide interactions."""

    # RNA-DNA hybrid nearest-neighbor parameters (37°C, 1M NaCl)
    # From Sugimoto et al. (1995) Biochemistry 34:11211-11216
    NN_PARAMS = {
        # Format: (ΔH kcal/mol, ΔS cal/mol·K)
        "rAA/dTT": (-7.8, -21.9),
        "rAC/dTG": (-5.9, -12.3),
        "rAG/dTC": (-9.1, -23.5),
        "rAU/dTA": (-8.3, -23.9),
        "rCA/dGT": (-9.0, -22.6),
        "rCC/dGG": (-9.3, -23.2),
        "rCG/dGC": (-16.3, -47.1),
        "rCU/dGA": (-7.0, -16.4),
        "rGA/dCT": (-5.5, -13.6),
        "rGC/dCG": (-8.0, -17.1),
        "rGG/dCC": (-12.8, -31.9),
        "rGU/dCA": (-7.8, -21.6),
        "rUA/dAT": (-7.8, -23.2),
        "rUC/dAG": (-8.6, -22.9),
        "rUG/dAC": (-10.4, -28.4),
        "rUU/dAA": (-11.5, -36.4),
    }

    # Mismatch penalty parameters (approximate)
    MISMATCH_PENALTY = {
        "transition": 2.5,  # A-G or C-T mismatch (kcal/mol)
        "transversion": 4.0,  # Other mismatches
        "gap": 5.0,  # Insertion/deletion
    }

    # Terminal penalties
    INIT_PENALTY = (0.0, -10.8)  # Initiation penalty (ΔH, ΔS)

    def __init__(
        self,
        na_conc: float = 0.05,  # Sodium concentration (M)
        mg_conc: float = 0.002,  # Magnesium concentration (M)
        oligo_conc: float = 1e-7,  # Oligo concentration (M)
        temperature: float = 37.0,
    ):  # Temperature (°C)
        """
        Initialize calculator with experimental conditions.

        Args:
            na_conc: Sodium ion concentration in M
            mg_conc: Magnesium ion concentration in M
            oligo_conc: Oligonucleotide concentration in M
            temperature: Temperature in Celsius
        """
        self.na_conc = na_conc
        self.mg_conc = mg_conc
        self.oligo_conc = oligo_conc
        self.temperature = temperature
        self.temp_kelvin = temperature + 273.15

        # Calculate effective salt concentration (Owczarzy et al. 2008)
        self.salt_correction = self._calculate_salt_correction()

    def _calculate_salt_correction(self) -> float:
        """Calculate salt correction factor for Tm."""
        # Monovalent salt correction
        na_factor = 12.5 * np.log10(self.na_conc)

        # Divalent salt correction (simplified)
        mg_factor = 3.0 * np.sqrt(self.mg_conc)

        return na_factor + mg_factor

    def _get_nn_params(self, rna_dimer: str) -> Tuple[float, float]:
        """
        Get nearest-neighbor parameters for RNA-DNA dimer.

        Args:
            rna_dimer: Two-base RNA sequence (e.g., 'AU')

        Returns:
            Tuple of (ΔH, ΔS) parameters
        """
        # Map RNA dimer to DNA complement
        complement = {"A": "T", "U": "A", "G": "C", "C": "G", "T": "A"}

        # Get DNA complement (reverse)
        dna_dimer = "".join(complement.get(b, "N") for b in rna_dimer[::-1])

        # Create lookup key
        key = f"r{rna_dimer}/d{dna_dimer}"

        # Return parameters or average if not found
        return self.NN_PARAMS.get(key, (-8.0, -22.0))

    def calculate_thermodynamics(
        self, guide_seq: str, target_seq: str, mismatches: Optional[List[int]] = None
    ) -> ThermodynamicProfile:
        """
        Calculate thermodynamic profile for guide RNA-target interaction.

        Args:
            guide_seq: Guide RNA sequence (5' to 3')
            target_seq: Target DNA sequence (5' to 3')
            mismatches: List of mismatch positions (0-indexed)

        Returns:
            ThermodynamicProfile with calculated parameters
        """
        guide_seq = guide_seq.upper().replace("T", "U")
        target_seq = target_seq.upper()

        # Initialize with initiation penalty
        total_dh, total_ds = self.INIT_PENALTY

        # Calculate nearest-neighbor contributions
        for i in range(len(guide_seq) - 1):
            dimer = guide_seq[i : i + 2]
            dh, ds = self._get_nn_params(dimer)
            total_dh += dh
            total_ds += ds

        # Add mismatch penalties
        mismatch_penalty = 0.0
        if mismatches:
            for pos in mismatches:
                if pos < len(guide_seq) and pos < len(target_seq):
                    g_base = guide_seq[pos]
                    t_base = target_seq[pos]

                    # Determine mismatch type
                    if (g_base in "AG" and t_base in "TC") or (
                        g_base in "UC" and t_base in "AG"
                    ):
                        mismatch_penalty += self.MISMATCH_PENALTY["transition"]
                    else:
                        mismatch_penalty += self.MISMATCH_PENALTY["transversion"]

        # Calculate ΔG at experimental temperature
        # ΔG = ΔH - TΔS
        delta_g = total_dh - (self.temp_kelvin * total_ds / 1000)
        delta_g += mismatch_penalty  # Add mismatch penalty to free energy

        # Calculate melting temperature
        # Tm = ΔH / (ΔS + R·ln(C))
        R = 1.987  # Gas constant (cal/mol·K)
        conc_term = R * np.log(self.oligo_conc / 4)
        tm_kelvin = (total_dh * 1000) / (total_ds + conc_term)
        tm_celsius = tm_kelvin - 273.15 + self.salt_correction

        # Calculate binding efficiency (simplified model)
        # Based on fraction bound at experimental temperature
        if delta_g < -30:
            binding_eff = 1.0
        elif delta_g > 0:
            binding_eff = 0.0
        else:
            # Sigmoid function normalized to ΔG range
            binding_eff = 1 / (1 + np.exp(delta_g / 5))

        return ThermodynamicProfile(
            sequence=guide_seq,
            target=target_seq,
            delta_g=delta_g,
            delta_h=total_dh,
            delta_s=total_ds,
            tm=tm_celsius,
            binding_efficiency=binding_eff,
            mismatch_penalty=mismatch_penalty,
        )

    def predict_off_target_binding(
        self, guide_seq: str, off_target_seq: str, pam_proximal_weight: float = 2.0
    ) -> float:
        """
        Predict off-target binding probability.

        Args:
            guide_seq: Guide RNA sequence
            off_target_seq: Potential off-target sequence
            pam_proximal_weight: Weight for PAM-proximal mismatches

        Returns:
            Off-target binding probability (0-1)
        """
        # Find mismatches
        mismatches = []
        for i, (g, t) in enumerate(zip(guide_seq, off_target_seq)):
            if g != t:
                mismatches.append(i)

        if not mismatches:
            return 1.0  # Perfect match

        # Calculate thermodynamic profile
        profile = self.calculate_thermodynamics(guide_seq, off_target_seq, mismatches)

        # Apply PAM-proximal weighting (last 10 bases are critical)
        pam_proximal_mismatches = sum(1 for m in mismatches if m >= len(guide_seq) - 10)

        # Reduce binding probability based on PAM-proximal mismatches
        if pam_proximal_mismatches > 0:
            profile.binding_efficiency *= 0.5 ** (
                pam_proximal_mismatches * pam_proximal_weight
            )

        return profile.binding_efficiency

    def calculate_specificity_score(
        self,
        guide_seq: str,
        on_target_profile: ThermodynamicProfile,
        off_target_profiles: List[ThermodynamicProfile],
    ) -> float:
        """
        Calculate specificity score based on thermodynamic profiles.

        Args:
            guide_seq: Guide RNA sequence
            on_target_profile: Thermodynamic profile for on-target
            off_target_profiles: List of off-target profiles

        Returns:
            Specificity score (0-1, higher is more specific)
        """
        if not off_target_profiles:
            return 1.0

        # Calculate ratio of on-target to off-target binding
        on_target_binding = on_target_profile.binding_efficiency

        # Sum of off-target binding probabilities
        off_target_sum = sum(p.binding_efficiency for p in off_target_profiles)

        # Specificity = on-target / (on-target + off-targets)
        specificity = on_target_binding / (on_target_binding + off_target_sum)

        return specificity


def analyze_guide_thermodynamics(
    guide_seq: str, target_seqs: Dict[str, str], conditions: Optional[Dict] = None
) -> Dict:
    """
    Comprehensive thermodynamic analysis of a guide RNA.

    Args:
        guide_seq: Guide RNA sequence
        target_seqs: Dictionary of target sequences {id: sequence}
        conditions: Experimental conditions (na_conc, mg_conc, etc.)

    Returns:
        Dictionary with thermodynamic analysis results
    """
    # Initialize calculator with conditions
    if conditions is None:
        conditions = {}

    calc = ThermodynamicCalculator(**conditions)

    results = {
        "guide_sequence": guide_seq,
        "profiles": {},
        "average_tm": 0,
        "average_delta_g": 0,
        "stability_scores": [],
        "binding_efficiencies": [],
    }

    # Analyze each target
    for target_id, target_seq in target_seqs.items():
        # Extract 20bp target region if longer
        if len(target_seq) > 20:
            # Assume target region is at the beginning for simplicity
            target_seq = target_seq[:20]

        profile = calc.calculate_thermodynamics(guide_seq, target_seq)
        results["profiles"][target_id] = profile

        results["stability_scores"].append(profile.stability_score)
        results["binding_efficiencies"].append(profile.binding_efficiency)

    # Calculate averages
    if results["profiles"]:
        results["average_tm"] = np.mean([p.tm for p in results["profiles"].values()])
        results["average_delta_g"] = np.mean(
            [p.delta_g for p in results["profiles"].values()]
        )
        results["mean_stability"] = np.mean(results["stability_scores"])
        results["mean_binding_efficiency"] = np.mean(results["binding_efficiencies"])

    return results
