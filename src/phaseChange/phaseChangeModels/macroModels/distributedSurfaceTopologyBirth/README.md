# distributedSurfaceTopologyBirth

Capability-first wall topology-birth model for resolved VOF boiling.

This model does **not** represent a calibrated cavity population. It supplies only
the minimum finite vapor topology required for an initially all-liquid VOF field
to begin resolved boiling. After a seed reaches its finite target or its creation
window expires, the model stops sourcing that location. Subsequent growth,
deformation, departure, motion, coalescence, evaporation, and condensation are
owned by the normal VOF/phase-change solver.

## Event process

Every configured heater-patch face is a potential event location. A face is
eligible when its adjacent cell is sufficiently liquid and its wall temperature
exceeds `Tsat + superheatThreshold`. The event hazard is

```
lambda = baseHazardRate
       * ((DeltaT - superheatThreshold)/hazardSuperheatScale)^hazardExponent
```

with optional `maxHazardRate` clipping. The per-step event probability is
`1 - exp(-lambda*Aface*deltaT)`, so the expected birth rate scales with physical
surface area and timestep rather than face count.

No permanent site coordinates are stored. A finite exclusion/cooldown footprint
around each event prevents immediate repeated insertion at the same location.
Existing vapor on the candidate wall-adjacent cell also suppresses a new event
through `minimumLiquidFraction`.

## Seed handoff

`seedRadiusCells` converts the local wall-adjacent cell volume into a physical
seed radius. `maximumSeedRadius` may cap that radius. The target is a smoothed
hemispherical vapor cap centered on the wall face. The alpha source is spread
over `seedCreationDuration`, capped by `maxAlphaVaporPerStep`, and accompanied
by liquid mass removal and a latent-energy sink. `thermalReserveFraction > 0`
optionally limits each cell's source by its instantaneous sensible superheat.

The source automatically falls to zero where the resolved vapor fraction has
already reached the finite target.

## Current qualification boundary

The first implementation is intentionally serial-only. Event ownership,
decomposition-independent sampling, state persistence across restart, and
cross-rank seed construction remain explicit MPI/restart work items.
