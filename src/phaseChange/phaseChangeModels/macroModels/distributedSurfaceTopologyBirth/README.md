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

If `exclusionRadius` is zero, the footprint scales with the created seed radius
through `exclusionRadiusFactor`. This is the preferred capability-mode setting
because it permits denser surface populations as the mesh/seed scale is reduced.

## MPI event synchronization

The MPI implementation preserves a physical whole-surface event process instead
of assigning permanent nucleation sites to processor domains.

Each rank evaluates only the heater faces it owns. The stateless stochastic draw
is keyed by time index, configured patch ordinal, and quantized physical face
centre, so it does not depend on processor-local face numbering. Candidate events
are gathered to rank 0, deterministically ordered by physical location, filtered
by the global `maxBirthsPerStep` and exclusion/cooldown history, then broadcast to
all ranks.

Every accepted physical event is therefore represented identically on all ranks.
Each rank builds only the portion of the hemispherical seed stencil that lies in
its local cells. The accepted-event history is replicated, so exclusion/cooldown
footprints cross processor boundaries. Source integrals and topology-defect
metrics are globally reduced before diagnostic output.

The first qualification target is decomposition invariance of event locations,
birth count, resolved alpha visibility, and integral diagnostics for 1, 4, 8,
and 16 ranks. Restart persistence of active event state is still a separate work
item.

## Seed handoff

`seedRadiusCells` converts the local wall-adjacent cell volume into a physical
seed radius. `maximumSeedRadius` may cap that radius. The target is a smoothed
hemispherical vapor cap centered on the wall face. The alpha source is spread
over `seedCreationDuration` and capped by `maxAlphaVaporPerStep`.

`sourceCouplingMode` controls how that finite alpha topology is coupled:

- `conservative`: retain the original alpha + liquid-mass + latent-energy source
  path. `thermalReserveFraction > 0` may limit the source by instantaneous
  sensible superheat. This mode remains available for later conservative-handoff
  qualification.
- `topologyOnly`: apply only the finite alpha topology. The model does **not** add
  its artificial mass source or latent sink. Instead it writes the equivalent
  mass/latent defect rate fields and accumulates the corresponding mass and
  energy defects in the `DSTB` diagnostics. Ordinary resolved phase change owns
  mass/energy transfer after the topology has been created.

The topology-only path is intentionally a quantified capability fallback, not a
claim of conservation. It exists to avoid the unphysical local energy impulse
observed when a mesh-resolved cap is forced into existence faster than the heater
can physically supply its latent energy.

The source automatically falls to zero where the resolved vapor fraction has
already reached the finite target.

## Diagnostic fields

With diagnostics enabled the model writes:

- `<prefix>Superheat`
- `<prefix>EligibleMask`
- `<prefix>ActiveSeedMask`
- `<prefix>AlphaSource`
- `<prefix>MassSource`
- `<prefix>LatentSink`
- `<prefix>EquivalentMassDefectRate`
- `<prefix>EquivalentLatentDefectRate`

In `topologyOnly` mode the actual `MassSource` and `LatentSink` fields remain
zero while the equivalent-defect fields quantify what the conservative seed
creation would have required.

## Current qualification boundary

MPI event synchronization and cross-rank seed construction are implemented but
must be qualified against the serial baseline before production use. In
particular, decomposition invariance, restart persistence, and longer multi-node
runs remain explicit verification items.
