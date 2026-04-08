.. _meld-algorithm:

==============
MELD Algorithm
==============

MELD (Modeling Employing Limited Data) is a Bayesian inference framework for
determining biomolecular structures from sparse, ambiguous, or noisy experimental
data [MacCallum2015]_, [Perez2015]_. This document describes how MELD works at a
conceptual level, focusing on the role of the ``alpha`` parameter and the scalers
that use it to modulate restraint strengths across a replica exchange ladder.

Overview
--------

MELD combines a physical force field with experimental restraints to sample
the posterior distribution

.. math::

    p(x | D) \propto p(x) \, p(D | x)

where :math:`x` is the molecular structure, :math:`p(x)` encodes the physical
prior (the force field), and :math:`p(D | x)` is the likelihood of the data
given the structure.

In practice, experimental data are often ambiguous or partially incorrect (e.g.,
NOE assignments that may be wrong, or predicted contacts that are not all
satisfied simultaneously). To accommodate this, MELD introduces *selectable
restraints* that allow only the most consistent subset of restraints to be active
at any step, and uses a *replica exchange* scheme to explore the posterior
efficiently.

The Alpha Parameter
-------------------

Each replica in a MELD simulation is assigned a value :math:`\alpha` that ranges
continuously from **0.0** (the lowest, most-physical replica) to **1.0** (the
highest, least-restrained replica):

- **alpha = 0**: restraints are applied at full strength. The simulation is most
  closely guided by the experimental data.
- **alpha = 1**: restraints are applied at (near) zero strength. The simulation
  runs essentially as a free, unbiased molecular dynamics simulation.

When there are :math:`N` replicas, the :math:`i`-th replica (zero-indexed) is
assigned

.. math::

    \alpha_i = \frac{i}{N - 1}, \quad i = 0, 1, \dots, N-1

The lowest replica (:math:`\alpha = 0`) samples the target posterior distribution.
The higher replicas act as "heat baths" that help the system escape local minima
and ensure global sampling. Replica exchange moves periodically swap configurations
between adjacent replicas, allowing structures that form at high :math:`\alpha` to
relax into the low-:math:`\alpha` regime where the full restraint forces are active.

The effective force constant for a restrained degree of freedom at a given replica
is

.. math::

    k_\mathrm{eff} = k \cdot S(\alpha)

where :math:`k` is the nominal (maximum) force constant and :math:`S(\alpha)` is
the *scaler* function described below.

Scalers
-------

Scalers are functions that map :math:`\alpha \in [0, 1]` to a scaling factor in
:math:`[0, 1]` (or more generally ``[strength_at_alpha_max,
strength_at_alpha_min]``). Every restraint in MELD may be assigned an optional
scaler; if none is provided, a :class:`~meld.system.scalers.ConstantScaler` is
used that always returns 1.0.

All scalers accept at minimum:

- ``alpha_min``: the value of :math:`\alpha` below which the scaling is at its
  maximum (``strength_at_alpha_min``).
- ``alpha_max``: the value of :math:`\alpha` above which the scaling is at its
  minimum (``strength_at_alpha_max``).
- ``strength_at_alpha_min``: strength (default 1.0) returned when
  :math:`\alpha \le` ``alpha_min``.
- ``strength_at_alpha_max``: strength (default ~0.001) returned when
  :math:`\alpha \ge` ``alpha_max``.

The following scalers are available.

ConstantScaler
~~~~~~~~~~~~~~

Always returns 1.0 regardless of :math:`\alpha`. Use this when a restraint should
be active at the same strength across all replicas.

.. code-block:: python

    scaler = system.restraints.create_scaler("constant")

LinearScaler
~~~~~~~~~~~~

Interpolates linearly between ``strength_at_alpha_min`` and
``strength_at_alpha_max`` across the range ``[alpha_min, alpha_max]``.

.. math::

    S(\alpha) = \begin{cases}
        s_\mathrm{min} & \alpha \le \alpha_\mathrm{min} \\[4pt]
        s_\mathrm{min} + (s_\mathrm{max} - s_\mathrm{min})
            \dfrac{\alpha - \alpha_\mathrm{min}}{\alpha_\mathrm{max} - \alpha_\mathrm{min}}
            & \alpha_\mathrm{min} < \alpha < \alpha_\mathrm{max} \\[4pt]
        s_\mathrm{max} & \alpha \ge \alpha_\mathrm{max}
    \end{cases}

where :math:`s_\mathrm{min}` = ``strength_at_alpha_min`` and
:math:`s_\mathrm{max}` = ``strength_at_alpha_max``.

.. code-block:: python

    scaler = system.restraints.create_scaler(
        "linear", alpha_min=0.0, alpha_max=0.5
    )

NonLinearScaler
~~~~~~~~~~~~~~~

Interpolates between the two strength values using an exponential (non-linear)
curve, giving finer control over how quickly restraints are turned off as
:math:`\alpha` increases. The ``factor`` parameter (``>= 1``) controls the
degree of non-linearity: larger values produce a more abrupt transition.

.. math::

    S(\alpha) = \frac{e^{f(1 - \delta)} - 1}{e^f - 1}

where :math:`\delta = (\alpha - \alpha_\mathrm{min}) /
(\alpha_\mathrm{max} - \alpha_\mathrm{min})` and :math:`f` = ``factor``.

.. code-block:: python

    scaler = system.restraints.create_scaler(
        "nonlinear", alpha_min=0.0, alpha_max=0.5, factor=4.0
    )

GeometricScaler
~~~~~~~~~~~~~~~

Interpolates between the two strength values on a geometric (logarithmic) scale.
This is useful when the strengths span several orders of magnitude and a uniform
distribution of the *ratio* of successive strengths is desired.

.. math::

    S(\alpha) = s_\mathrm{min} \exp\!\left(
        \frac{\alpha - \alpha_\mathrm{min}}{\alpha_\mathrm{max} - \alpha_\mathrm{min}}
        \ln\frac{s_\mathrm{max}}{s_\mathrm{min}}
    \right)

.. code-block:: python

    scaler = system.restraints.create_scaler(
        "geometric", alpha_min=0.0, alpha_max=1.0,
        strength_at_alpha_min=1.0, strength_at_alpha_max=1e-3
    )

Plateau Scalers
~~~~~~~~~~~~~~~

Three plateau-shaped scalers are available. They all ramp up from
``strength_at_alpha_max`` to ``strength_at_alpha_min`` between ``alpha_min`` and
``alpha_one``, hold at ``strength_at_alpha_min`` between ``alpha_one`` and
``alpha_two``, and then ramp back down to ``strength_at_alpha_max`` between
``alpha_two`` and ``alpha_max``. This shape is useful when a restraint should be
active only in a middle band of replicas.

The three variants differ in how they interpolate in the ramp regions:

- **PlateauLinearScaler** (``"plateau"``): linear ramps.
- **PlateauNonLinearScaler** (``"plateaunonlinear"``): exponential ramps,
  controlled by a ``factor`` parameter.
- **PlateauSmoothScaler** (``"plateausmooth"``): smooth (cubic Hermite)
  ramps with zero derivative at the transition points.

.. code-block:: python

    scaler = system.restraints.create_scaler(
        "plateau",
        alpha_min=0.0, alpha_one=0.3, alpha_two=0.7, alpha_max=1.0
    )

Ramps
-----

Ramps are analogous to scalers but map *simulation time* (step number) to a
scaling factor rather than :math:`\alpha`. They are typically used to slowly turn
on restraint forces at the beginning of a simulation to avoid instabilities from
large initial forces. Available ramps are:

- **ConstantRamp** (``"constant_ramp"``): always returns 1.0.
- **LinearRamp** (``"linear_ramp"``): linearly interpolates between
  ``start_weight`` and ``end_weight`` over a specified time window.
- **NonLinearRamp** (``"nonlinear_ramp"``): interpolates non-linearly, with
  a ``factor`` controlling the shape of the curve.
- **TimeRampSwitcher** (``"ramp_switcher"``): uses one ramp before a
  ``switching_time`` and a second ramp afterwards.

Ramps are created with the same :meth:`~meld.system.restraints.RestraintManager.create_scaler`
method as scalers:

.. code-block:: python

    ramp = system.restraints.create_scaler(
        "linear_ramp", start_time=0, end_time=200,
        start_weight=0.0, end_weight=1.0
    )
    r = system.restraints.create_restraint(rest_key, ramp=ramp, ...)

Positioners
-----------

Positioners are a third type of alpha-dependent mapper. Instead of returning a
dimensionless scale factor in ``[0, 1]``, they return a *position* (e.g., a
distance in nanometers) within a defined range. This allows the equilibrium values
of restraints to change smoothly across the replica ladder.

- **ConstantPositioner** (``"constant_positioner"``): always returns the same
  value.
- **LinearPositioner** (``"linear_positioner"``): linearly interpolates between
  ``pos_min`` and ``pos_max`` as :math:`\alpha` goes from ``alpha_min`` to
  ``alpha_max``.

.. code-block:: python

    positioner = system.restraints.create_scaler(
        "linear_positioner",
        alpha_min=0.0, alpha_max=1.0,
        pos_min=0.3 * u.nanometer,
        pos_max=1.0 * u.nanometer,
    )
    r = system.restraints.create_restraint(
        "distance", r3=positioner, ...
    )

Combining Scalers, Ramps, and Positioners
-----------------------------------------

A single restraint may have:

1. A **scaler** — controls the force constant as a function of :math:`\alpha`.
2. A **ramp** — modulates the force constant during the initial equilibration.
3. **Positioners** — replace fixed distance (or angle) parameters with
   :math:`\alpha`-dependent values.

The effective force constant at step :math:`t` in replica :math:`i` is

.. math::

    k_\mathrm{eff} = k \cdot S(\alpha_i) \cdot R(t)

where :math:`S` is the scaler and :math:`R` is the ramp.

Example
-------

The following example creates a distance restraint that is fully active at low
:math:`\alpha` replicas and nearly off at high :math:`\alpha` replicas, and is
gently ramped on over the first 200 steps of the simulation:

.. code-block:: python

    scaler = system.restraints.create_scaler(
        "nonlinear", alpha_min=0.0, alpha_max=0.5, factor=4.0
    )
    ramp = system.restraints.create_scaler(
        "linear_ramp", start_time=0, end_time=200,
        start_weight=0.0, end_weight=1.0
    )
    r = system.restraints.create_restraint(
        "distance",
        scaler=scaler,
        ramp=ramp,
        atom1=atom_i,
        atom2=atom_j,
        r1=0.0 * u.nanometer,
        r2=0.2 * u.nanometer,
        r3=0.5 * u.nanometer,
        r4=0.7 * u.nanometer,
        k=250.0 * u.kilojoule_per_mole / u.nanometer**2,
    )
    system.restraints.add_as_always_active(r)

References
----------

.. [MacCallum2015] J.L. MacCallum, A. Perez, and K.A. Dill, Determining protein
   structures by combining semireliable data with atomistic physical models by
   Bayesian inference, *PNAS*, 2015, **112** (22), pp. 6985–6990.

.. [Perez2015] A. Perez, J.L. MacCallum, and K.A. Dill, Accelerating molecular
   simulations of proteins using Bayesian inference on weak information, *PNAS*,
   2015, **112** (38), pp. 11846–11851.
