// Case224: field and boundary diagnostics. Collect one reference inventory
// BEFORE the first time step; never advance the reference in a PIMPLE loop.
const Switch massBalanceDiagnostics
(
    runTime.controlDict().lookupOrDefault<Switch>
    (
        "massBalanceDiagnostics",
        false
    )
);
const label massBalanceReportInterval =
    runTime.controlDict().lookupOrDefault<label>
    (
        "massBalanceReportInterval",
        1
    );
if (massBalanceDiagnostics && massBalanceReportInterval < 1)
{
    FatalErrorInFunction
        << "massBalanceReportInterval must be >= 1"
        << exit(FatalError);
}

scalarField massBalanceInitialMass(fluidRegions.size(), 0.0);
scalarField massBalanceLastMass(fluidRegions.size(), 0.0);
scalarField massBalanceCumulativeAdvOut(fluidRegions.size(), 0.0);
scalarField massBalanceCumulativeCorrectedOut(fluidRegions.size(), 0.0);
scalarField massBalanceCumulativeTopAdvOut(fluidRegions.size(), 0.0);
scalarField massBalanceCumulativeTopCorrectedOut(fluidRegions.size(), 0.0);
scalarField massBalanceCumulativeMacroBirthMass(fluidRegions.size(), 0.0);
scalarField massBalanceCumulativeMacroBirthVolume(fluidRegions.size(), 0.0);
label massBalanceStepCount = 0;

if (massBalanceDiagnostics)
{
    forAll(fluidRegions, mbRegionI)
    {
        const scalar initialInventory =
            fvc::domainIntegrate(rhoFluid[mbRegionI]).value();
        const scalar initialVaporVolume =
            fvc::domainIntegrate(thermoFluid[mbRegionI].alpha2()).value();
        massBalanceInitialMass[mbRegionI] = initialInventory;
        massBalanceLastMass[mbRegionI] = initialInventory;
        if (Pstream::master())
        {
            Info<< "MASSCHK_INIT"
                << " region=" << fluidRegions[mbRegionI].name()
                << " time=" << runTime.value()
                << " massKg=" << initialInventory
                << " vaporVolumeM3=" << initialVaporVolume
                << endl;
        }
    }
}
