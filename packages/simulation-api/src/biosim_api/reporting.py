# SAVE TO: packages\simulation-api\src\biosim_api\reporting.py
"""
Expert-metric and lay-summary computation for SimulationRunner.

Moved verbatim out of runner.py (2026-10-20) to keep that module readable;
SimulationRunner inherits these methods, so behaviour is unchanged.
"""

from __future__ import annotations

import numpy as np

from biosim_ode import (
    PathwayParameters,
)

from .output import ExpertMetrics, LaySummary

# Threshold on ConfidenceBand("G").mean_relative_width() above which a
# Monte Carlo UQ ensemble's glucose trajectory is considered "wide enough"
# to flag in the lay summary, even when the point-estimate Boolean watchdog
# reports high confidence. An engineering choice (same category as
# MIN_GB_THRESHOLD_MARGIN_MM/ConfidenceBand.mean_relative_width()'s own
# floor_fraction), not a literature-derived statistic -- physiologically
# real species have historically measured ~0.03-0.08 under normal
# conditions in this platform's own diagnostics (see result.py's docstring
# history), so 0.15 is comfortably above ordinary noise while still
# catching a genuinely destabilized ensemble.
UQ_WIDE_GLUCOSE_RELATIVE_WIDTH = 0.15


class ReportingMixin:
    """
    Turns a finished simulation into ExpertMetrics and a LaySummary.
    """

    # ------------------------------------------------------------------
    # Expert metrics computation
    # ------------------------------------------------------------------

    def _compute_expert_metrics(self, sim, params: PathwayParameters) -> ExpertMetrics:
        if not sim.success or not sim.species:
            return ExpertMetrics(
                glucose_peak_uM=0, glucose_peak_time_s=0, glucose_auc_uM_s=0,
                glucose_baseline_uM=params.Gb * 1000.0,
                glucose_nadir_uM=params.Gb * 1000.0,
                insulin_peak_uM=0, insulin_peak_time_s=0,
                mtor_peak=0, ampk_nadir=1.0,
                time_to_glucose_return_s=None,
                overall_watchdog_confidence=0.0,
                n_ode_evaluations=0,
            )

        G = sim.species["G"]
        I = sim.species["I"]

        peak_G, peak_G_t = sim.glucose_peak()
        peak_I, peak_I_t = sim.insulin_peak()
        baseline_G = float(G[0])
        nadir_G = float(np.min(G))
        auc_G = float(np.trapezoid(G, sim.times_s))
        mtor_peak = float(np.max(sim.species["mTOR"]))
        ampk_nadir = float(np.min(sim.species["AMPK"]))

        # iAUC(0-5h): the standard clinical incremental-AUC metric (see
        # ExpertMetrics.glucose_iauc_5h_uM_s docstring). 0.0 if the sim
        # didn't run long enough for a meaningful 5h window.
        FIVE_HOURS_S = 18000.0
        iauc_5h = 0.0
        if sim.times_s[-1] >= FIVE_HOURS_S:
            mask_5h = sim.times_s <= FIVE_HOURS_S
            incremental = np.clip(G[mask_5h] - baseline_G, 0.0, None)
            iauc_5h = float(np.trapezoid(incremental, sim.times_s[mask_5h]))

        # iAUC(0-4h): same metric, fixed to a 4h window -- the standard
        # "mixed meal test" endpoint duration (e.g. Hansen et al. 2024,
        # PMID:39049802, for metformin). 0.0 if the sim didn't run long
        # enough for a meaningful 4h window.
        FOUR_HOURS_S = 14400.0
        iauc_4h = 0.0
        if sim.times_s[-1] >= FOUR_HOURS_S:
            mask_4h = sim.times_s <= FOUR_HOURS_S
            incremental_4h = np.clip(G[mask_4h] - baseline_G, 0.0, None)
            iauc_4h = float(np.trapezoid(incremental_4h, sim.times_s[mask_4h]))

        return_t = sim.time_to_return_to_baseline("G", baseline_G, threshold_fraction=0.10)

        confidence = 1.0
        if sim.watchdog:
            from biosim_ode.system import ODESystem as _ODE
            ode = _ODE()
            confidence = ode.overall_confidence(sim)

        return ExpertMetrics(
            glucose_peak_uM=peak_G,
            glucose_peak_time_s=peak_G_t,
            glucose_auc_uM_s=auc_G,
            glucose_baseline_uM=baseline_G,
            glucose_nadir_uM=nadir_G,
            insulin_peak_uM=peak_I,
            insulin_peak_time_s=peak_I_t,
            mtor_peak=mtor_peak,
            ampk_nadir=ampk_nadir,
            time_to_glucose_return_s=return_t,
            overall_watchdog_confidence=confidence,
            n_ode_evaluations=sim.solver_nfev,
            glucose_iauc_5h_uM_s=iauc_5h,
            glucose_iauc_4h_uM_s=iauc_4h,
        )

    # ------------------------------------------------------------------
    # Lay summary computation
    # ------------------------------------------------------------------

    def _compute_lay_summary(self, sim, metrics: ExpertMetrics, uq) -> LaySummary:
        if not sim.success:
            empty = np.zeros(10)
            return LaySummary(
                times_min=np.linspace(0, 240, 10),
                energy_level=empty + 5.0,
                hunger_level=empty + 5.0,
                energy_peak_min=0.0,
                energy_dip_min=None,
                has_energy_dip=False,
                glucose_narrative="Simulation did not complete successfully.",
                energy_narrative="",
                confidence_note="Results unavailable.",
                is_high_confidence=False,
            )

        G = sim.species["G"]
        I = sim.species["I"]
        times_min = sim.times_s / 60.0
        G_baseline = float(G[0])

        # Energy level proxy (0-10)
        # Rises with glucose above baseline, drops when below
        G_norm = (G - G_baseline) / G_baseline   # fractional deviation
        energy = 5.0 + 3.0 * G_norm - 2.0 * np.maximum(-G_norm, 0)
        energy = np.clip(energy, 0.0, 10.0)

        # Hunger proxy (0-10)
        hunger = 5.0 - 3.0 * G_norm + 2.0 * np.maximum(-G_norm, 0)
        hunger = np.clip(hunger, 0.0, 10.0)

        # Energy peak timing
        energy_peak_idx = int(np.argmax(energy))
        energy_peak_min = float(times_min[energy_peak_idx])

        # Energy dip: look for a drop below starting level after the peak
        post_peak_energy = energy[energy_peak_idx:]
        dip_idx = None
        if len(post_peak_energy) > 0:
            min_post_peak = float(np.min(post_peak_energy))
            if min_post_peak < 4.5:  # meaningful dip below baseline
                dip_offset = int(np.argmin(post_peak_energy))
                dip_idx = energy_peak_idx + dip_offset

        has_dip = dip_idx is not None
        dip_min = float(times_min[dip_idx]) if dip_idx is not None else None

        # Narrative generation
        peak_mM = round(metrics.glucose_peak_mM, 1)
        peak_min = round(metrics.glucose_peak_time_s / 60.0, 0)
        glucose_narrative = (
            f"Your blood sugar is estimated to peak at around {peak_mM} mmol/L "
            f"approximately {peak_min:.0f} minutes after eating."
        )
        if metrics.time_to_glucose_return_s:
            return_min = round(metrics.time_to_glucose_return_s / 60.0, 0)
            baseline_mM = round(metrics.glucose_baseline_uM / 1000.0, 1)
            # Regression: time_to_glucose_return_s measures return to this
            # patient's OWN baseline (glucose_baseline_uM = G[0]), not a
            # population-normal value -- unconditionally calling that
            # "toward normal" is false whenever the baseline itself is
            # elevated (e.g. the T2D disease-state preset's Gb x1.3, ~6.5
            # mM, already in the ADA 2023 pre-diabetic fasting range of
            # 5.6-6.9 mM). Same fix already applied to
            # biosim_lay.report.py's sibling narrative for the identical
            # metric -- this is a separate code path (this package does
            # not, and should not, depend on biosim_lay) that the earlier
            # fix did not reach, despite being exposed through the same
            # REST API response (LaySummaryOut.glucose_narrative).
            if baseline_mM >= 5.6:
                glucose_narrative += (
                    f" It is expected to return to your own baseline ({baseline_mM} mmol/L) "
                    f"within {return_min:.0f} minutes. Note: this baseline itself is above "
                    f"the normal fasting range (<5.6 mmol/L)."
                )
            else:
                glucose_narrative += f" It is expected to return toward normal within {return_min:.0f} minutes."
        elif metrics.time_to_glucose_return_s is None and len(times_min) > 0:
            # Not back within the simulated window: say so rather than stop
            # after the peak (reachable for a T2D 75 g meal since B7).
            baseline_mM = round(metrics.glucose_baseline_uM / 1000.0, 1)
            window_h = times_min[-1] / 60.0
            glucose_narrative += (
                f" It has not returned to your own baseline ({baseline_mM} mmol/L) "
                f"by the end of the {window_h:g} {'hour' if window_h == 1 else 'hours'} simulated."
            )

        if has_dip:
            energy_narrative = (
                f"You may notice a rise in energy around {energy_peak_min:.0f} minutes, "
                f"followed by an energy dip around {dip_min:.0f} minutes as blood sugar "
                f"comes back down. This is sometimes called a 'sugar crash'."
            )
        else:
            energy_narrative = (
                f"Energy is expected to rise gradually, peaking around "
                f"{energy_peak_min:.0f} minutes, then slowly returning to normal "
                f"without a significant dip."
            )

        is_high_conf = metrics.overall_watchdog_confidence >= 0.8

        # Regression: is_high_confidence/confidence_note used to be derived
        # solely from the point-estimate Boolean watchdog score, structurally
        # blind to whether a Monte Carlo UQ ensemble (when the caller opted
        # in via ScenarioDefinition.run_uncertainty) actually produced wide
        # confidence bands or a high sample-failure rate. `uq` has been an
        # accepted-but-unused parameter here since before this ODE engine's
        # entire uncertainty-package overhaul -- every UQ correctness fix
        # (pathway_factory composition, the mean_relative_width() near-zero
        # floor, etc.) improved the ensemble's own *content* without this
        # lay-facing signal ever looking at it. A scenario could report
        # "high confidence" from a clean point estimate while its own UQ
        # ensemble shows the parameter uncertainty genuinely destabilizes
        # the glucose trajectory, or largely fails to solve at all.
        uq_note = ""
        if uq is not None:
            if not uq.is_reliable:
                is_high_conf = False
                uq_note = (
                    f" The uncertainty-quantification ensemble also had "
                    f"reduced reliability ({uq.success_rate:.0%} of parameter "
                    f"samples solved successfully, overall ensemble confidence "
                    f"{uq.overall_confidence:.2f})."
                )
            elif "G" in uq.bands and uq.bands["G"].mean_relative_width() > UQ_WIDE_GLUCOSE_RELATIVE_WIDTH:
                is_high_conf = False
                uq_note = (
                    " The uncertainty-quantification ensemble also shows a "
                    "wide range of plausible glucose trajectories for this "
                    "scenario under parameter uncertainty -- treat the exact "
                    "numbers above as illustrative, not precise."
                )

        confidence_note = (
            "This simulation has high confidence."
            if is_high_conf
            else "This simulation has moderate uncertainty. Treat results as approximate."
        ) + uq_note

        return LaySummary(
            times_min=times_min,
            energy_level=energy,
            hunger_level=hunger,
            energy_peak_min=energy_peak_min,
            energy_dip_min=dip_min,
            has_energy_dip=has_dip,
            glucose_narrative=glucose_narrative,
            energy_narrative=energy_narrative,
            confidence_note=confidence_note,
            is_high_confidence=is_high_conf,
        )
