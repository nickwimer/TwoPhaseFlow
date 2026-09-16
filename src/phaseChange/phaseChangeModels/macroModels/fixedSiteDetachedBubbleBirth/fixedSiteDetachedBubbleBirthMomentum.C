/*---------------------------------------------------------------------------*\
    Optional coarse-grid departure momentum closure for fixed-site boiling
\*---------------------------------------------------------------------------*/

#include "fixedSiteDetachedBubbleBirth.H"
#include "fvCFD.H"
#include "PstreamReduceOps.H"

void Foam::fixedSiteDetachedBubbleBirth::momentumSource
(
    fvVectorMatrix& UEqn
)
{
    const scalar targetUz = modelDict().lookupOrDefault<scalar>
    (
        "departureLaunchVelocity",
        0.0
    );

    // Zero is the backward-compatible default used by Case218 and all earlier
    // fixed-site qualification cases.
    if (targetUz <= SMALL)
    {
        if (targetUz < -SMALL)
        {
            FatalErrorInFunction
                << "departureLaunchVelocity must be >= 0; found "
                << targetUz << exit(FatalError);
        }
        return;
    }

    const scalar launchRelaxation = modelDict().lookupOrDefault<scalar>
    (
        "departureLaunchRelaxation",
        5.0e-5
    );
    if (launchRelaxation <= SMALL)
    {
        FatalErrorInFunction
            << "departureLaunchRelaxation must be > 0 when the launch closure "
            << "is enabled; found " << launchRelaxation
            << exit(FatalError);
    }

    // Assemble the same finite handoff state used by alpha/mass/energy source
    // terms. updateSources() is guarded by timeIndex, so this does not advance
    // a site twice when other phase-change hooks are evaluated later.
    updateSources();

    const fvMesh& mesh = phase1_.mesh();
    const scalar dt = max(mesh.time().deltaTValue(), scalar(VSMALL));
    const scalar relaxationTime = max(launchRelaxation, dt);
    const volVectorField& U = UEqn.psi();

    tmp<volScalarField> tRhoL = phase1_.thermo().rho();
    const volScalarField& rhoL = tRhoL();
    tmp<volScalarField> tRhoV = phase2_.thermo().rho();
    const volScalarField& rhoV = tRhoV();

    vectorField& source = UEqn.source();

    label localActiveCells = 0;
    scalar localIntegratedForceZ = 0;
    scalar localMaxWeight = 0;
    scalar localMinUz = GREAT;
    scalar localMaxUz = -GREAT;

    forAll(activeReleaseMask_, celli)
    {
        const scalar activeWeight = min
        (
            max(activeReleaseMask_[celli], scalar(0)),
            scalar(1)
        );
        if (activeWeight <= SMALL)
        {
            continue;
        }

        const scalar alphaLiquid = min
        (
            max(phase1_[celli], scalar(0)),
            scalar(1)
        );
        const scalar alphaVapor = scalar(1) - alphaLiquid;

        // The momentum closure follows the vapor being handed to VOF instead
        // of accelerating the entire stencil as liquid. During the first
        // source update alphaBirthSource*dt represents the scheduled new vapor
        // fraction before the alpha equation has applied it; on later updates
        // the resolved vapor fraction naturally becomes the larger weight.
        const scalar scheduledVapor = max
        (
            alphaBirthSource_[celli]*dt,
            scalar(0)
        );
        const scalar vaporWeight = min
        (
            max(alphaVapor, scheduledVapor),
            scalar(1)
        );
        const scalar weight = activeWeight*vaporWeight;
        if (weight <= SMALL)
        {
            continue;
        }

        const scalar deltaUz = max
        (
            targetUz - U[celli].z(),
            scalar(0)
        );
        if (deltaUz <= SMALL)
        {
            continue;
        }

        const scalar rhoMix =
            alphaLiquid*rhoL[celli]
          + (scalar(1) - alphaLiquid)*rhoV[celli];
        const scalar forceZ =
            rhoMix
           *mesh.V()[celli]
           *weight
           *deltaUz/relaxationTime;

        source[celli] += vector(0, 0, forceZ);
        localIntegratedForceZ += forceZ;
        localMaxWeight = max(localMaxWeight, weight);
        localMinUz = min(localMinUz, U[celli].z());
        localMaxUz = max(localMaxUz, U[celli].z());
        ++localActiveCells;
    }

    reduce(localActiveCells, sumOp<label>());
    reduce(localIntegratedForceZ, sumOp<scalar>());
    reduce(localMaxWeight, maxOp<scalar>());
    reduce(localMinUz, minOp<scalar>());
    reduce(localMaxUz, maxOp<scalar>());

    if (localActiveCells > 0 && Pstream::master())
    {
        Info<< "FIXED_SITE_LAUNCH"
            << " time=" << mesh.time().value()
            << " activeCells=" << localActiveCells
            << " targetUz=" << targetUz
            << " relaxationTime=" << relaxationTime
            << " integratedForceZ=" << localIntegratedForceZ
            << " maxWeight=" << localMaxWeight
            << " minUz=" << localMinUz
            << " maxUz=" << localMaxUz
            << endl;
    }
}

// ************************************************************************* //
