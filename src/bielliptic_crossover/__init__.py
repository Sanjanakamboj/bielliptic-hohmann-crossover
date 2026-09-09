"""Bi-elliptic vs. Hohmann transfer crossover analysis.

Release 1.0: the transfer equations and solvers (M2) plus the trade-study layer in
:mod:`bielliptic_crossover.trade` (crossover map, break-even locus, delta-v/time
trade, scoped recommendation). The analytical development and the verification
plan they implement are in ``DESIGN.md``.

Normalization::

    R      = r2 / r1            final-to-initial RADIUS ratio,      R > 1
    B      = rb / r1            intermediate apoapsis RADIUS ratio, B > R
    v1     = sqrt(mu / r1)      initial circular speed
    dv_bar = dv / v1            normalized delta-v

Quick start::

    >>> from bielliptic_crossover import threshold_R1, threshold_R2, classify_radius_ratio
    >>> round(threshold_R1(), 10)
    11.9387654726
    >>> round(threshold_R2(), 10)
    15.5817187388
    >>> classify_radius_ratio(20.0)
    'all_bielliptic'

No crossover constant is hardcoded anywhere: both thresholds are recomputed from
the transfer equations and cross-checked against exact polynomial forms.

The trade layer is imported explicitly as ``bielliptic_crossover.trade`` so that
this module's M2 surface stays exactly as verified.

Scope: idealized coplanar two-body impulsive transfers. No perturbations, finite
burns, plane changes or mission optimization -- see ``DESIGN.md`` for the full
limitations list.
"""

from __future__ import annotations

from .bielliptic import (
    bielliptic_burns_normalized,
    bielliptic_infinite_limit_normalized,
    bielliptic_total_derivative_wrt_B,
    bielliptic_total_normalized,
    burn_directions,
)
from .constants import (
    H1_REFERENCE,
    MU_EARTH,
    R1_REFERENCE,
    R_EARTH,
    V1_REFERENCE,
    altitude_from_radius,
    radius_from_altitude,
)
from .crossover import (
    B_RATIO_RESOLVABLE_MAX,
    REGIME_ALL_BIELLIPTIC,
    REGIME_HOHMANN_ONLY,
    REGIME_LARGE_B_BIELLIPTIC,
    BreakEvenResult,
    BStructure,
    bielliptic_infimum_normalized,
    break_even_B,
    classify_radius_ratio,
    scan_B_structure,
    threshold_R1,
    threshold_R1_from_polynomial,
    threshold_R2,
    threshold_R2_from_cubic,
)
from .dimensional import (
    bielliptic_burns_direct,
    bielliptic_total_dv,
    circular_speed,
    hohmann_burns_direct,
    hohmann_total_dv_direct,
    vis_viva_speed,
)
from .hohmann import (
    hohmann_burns_normalized,
    hohmann_total_dv,
    hohmann_total_normalized,
)
from .timing import (
    bielliptic_time_excess_at_B_equals_R,
    bielliptic_transfer_time,
    bielliptic_transfer_time_normalized,
    hohmann_transfer_time,
    hohmann_transfer_time_normalized,
    time_scale,
)

__version__ = "1.0.0"

__all__ = [
    "__version__",
    # constants
    "MU_EARTH",
    "R_EARTH",
    "H1_REFERENCE",
    "R1_REFERENCE",
    "V1_REFERENCE",
    "altitude_from_radius",
    "radius_from_altitude",
    # hohmann
    "hohmann_burns_normalized",
    "hohmann_total_normalized",
    "hohmann_total_dv",
    # bi-elliptic
    "bielliptic_burns_normalized",
    "bielliptic_total_normalized",
    "bielliptic_infinite_limit_normalized",
    "bielliptic_total_derivative_wrt_B",
    "burn_directions",
    # dimensional (independent vis-viva path)
    "vis_viva_speed",
    "circular_speed",
    "hohmann_burns_direct",
    "hohmann_total_dv_direct",
    "bielliptic_burns_direct",
    "bielliptic_total_dv",
    # timing
    "time_scale",
    "hohmann_transfer_time_normalized",
    "bielliptic_transfer_time_normalized",
    "hohmann_transfer_time",
    "bielliptic_transfer_time",
    "bielliptic_time_excess_at_B_equals_R",
    # crossover
    "threshold_R1",
    "threshold_R1_from_polynomial",
    "threshold_R2",
    "threshold_R2_from_cubic",
    "break_even_B",
    "BreakEvenResult",
    "classify_radius_ratio",
    "scan_B_structure",
    "BStructure",
    "B_RATIO_RESOLVABLE_MAX",
    "bielliptic_infimum_normalized",
    "REGIME_HOHMANN_ONLY",
    "REGIME_LARGE_B_BIELLIPTIC",
    "REGIME_ALL_BIELLIPTIC",
]
