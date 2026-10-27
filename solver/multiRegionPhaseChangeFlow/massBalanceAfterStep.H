// Case224: one sample per completed full time step (AFTER all outer
// correctors and all regions). Positive boundary flux is outward.
// IMPORTANT: the alpha-transport rhoPhi is evaluated before the final
// pressure correction, while phi is pressure-corrected later. These are two
// distinct estimators and NEITHER is guaranteed to be the solver's exact
// conservative total-mixture mass flux.
++massBalanceStepCount;
const scalar massBalanceDt = runTime.deltaTValue();
forAll(fluidRegions, mbRegionI)
{
    const fvMesh& mbMesh = fluidRegions[mbRegionI];
    const volScalarField& mbRho = rhoFluid[mbRegionI];
    const volScalarField& mbAlphaV = thermoFluid[mbRegionI].alpha2();
    const volScalarField& mbRhoV =
        thermoFluid[mbRegionI].thermo2().rho();
    const volScalarField& mbPhaseMassSource =
        PhaseChangeModelFluid[mbRegionI].massSource();
    const volScalarField& mbResolvedPsi0 =
        PhaseChangeModelFluid[mbRegionI].psi0();

    const word mbMacroMassName("populationSiteReleaseMassSource");
    const word mbMacroAlphaName("populationSiteReleaseAlphaSource");
    const volScalarField* mbMacroMassPtr = nullptr;
    const volScalarField* mbMacroAlphaPtr = nullptr;
    if (mbMesh.foundObject<volScalarField>(mbMacroMassName))
    {
        mbMacroMassPtr =
            &mbMesh.lookupObject<volScalarField>(mbMacroMassName);
    }
    if (mbMesh.foundObject<volScalarField>(mbMacroAlphaName))
    {
        mbMacroAlphaPtr =
            &mbMesh.lookupObject<volScalarField>(mbMacroAlphaName);
    }

    scalar mbLocalMass = 0.0;
    scalar mbLocalVaporVolume = 0.0;
    scalar mbLocalVaporMass = 0.0;
    scalar mbLocalPhaseSource = 0.0;
    scalar mbLocalResolvedPsi0 = 0.0;
    scalar mbLocalMacroBirthMass = 0.0;
    scalar mbLocalMacroBirthVolume = 0.0;
    const scalarField& mbCellVolumes = mbMesh.V();
    forAll(mbCellVolumes, mbCellI)
    {
        const scalar volume = mbCellVolumes[mbCellI];
        mbLocalMass += mbRho[mbCellI]*volume;
        mbLocalVaporVolume += mbAlphaV[mbCellI]*volume;
        mbLocalVaporMass +=
            mbAlphaV[mbCellI]*mbRhoV[mbCellI]*volume;
        mbLocalPhaseSource += mbPhaseMassSource[mbCellI]*volume;
        mbLocalResolvedPsi0 += mbResolvedPsi0[mbCellI]*volume;
        if (mbMacroMassPtr)
        {
            mbLocalMacroBirthMass +=
                (*mbMacroMassPtr)[mbCellI]*volume;
        }
        if (mbMacroAlphaPtr)
        {
            mbLocalMacroBirthVolume +=
                (*mbMacroAlphaPtr)[mbCellI]*volume;
        }
    }

    const surfaceScalarField& mbRhoPhi = rhoPhiFluid[mbRegionI];
    const surfaceScalarField& mbPhi = phiFluid[mbRegionI];
    tmp<surfaceScalarField> mbInterpolatedRho = fvc::interpolate(mbRho);
    const surfaceScalarField& mbFaceRho = mbInterpolatedRho();
    scalar mbLocalAdvOut = 0.0;
    scalar mbLocalCorrectedOut = 0.0;
    scalar mbLocalTopAdvOut = 0.0;
    scalar mbLocalTopCorrectedOut = 0.0;

    forAll(mbMesh.boundary(), mbPatchI)
    {
        const fvPatch& mbPatch = mbMesh.boundary()[mbPatchI];
        // Exclude processor and cyclic/coupled interfaces from EXTERNAL
        // fluxes. Empty physical patches have zero faces.
        if (mbPatch.coupled())
        {
            continue;
        }
        const auto& mbAdvPatch = mbRhoPhi.boundaryField()[mbPatchI];
        const auto& mbPhiPatch = mbPhi.boundaryField()[mbPatchI];
        const auto& mbRhoPatch = mbFaceRho.boundaryField()[mbPatchI];
        forAll(mbPhiPatch, mbFaceI)
        {
            const scalar adv = mbAdvPatch[mbFaceI];
            const scalar corrected =
                mbPhiPatch[mbFaceI]*mbRhoPatch[mbFaceI];
            mbLocalAdvOut += adv;
            mbLocalCorrectedOut += corrected;
            if (mbPatch.name() == "poolTop")
            {
                mbLocalTopAdvOut += adv;
                mbLocalTopCorrectedOut += corrected;
            }
        }
    }

    reduce(mbLocalMass, sumOp<scalar>());
    reduce(mbLocalVaporVolume, sumOp<scalar>());
    reduce(mbLocalVaporMass, sumOp<scalar>());
    reduce(mbLocalPhaseSource, sumOp<scalar>());
    reduce(mbLocalResolvedPsi0, sumOp<scalar>());
    reduce(mbLocalMacroBirthMass, sumOp<scalar>());
    reduce(mbLocalMacroBirthVolume, sumOp<scalar>());
    reduce(mbLocalAdvOut, sumOp<scalar>());
    reduce(mbLocalCorrectedOut, sumOp<scalar>());
    reduce(mbLocalTopAdvOut, sumOp<scalar>());
    reduce(mbLocalTopCorrectedOut, sumOp<scalar>());

    massBalanceCumulativeAdvOut[mbRegionI] +=
        massBalanceDt*mbLocalAdvOut;
    massBalanceCumulativeCorrectedOut[mbRegionI] +=
        massBalanceDt*mbLocalCorrectedOut;
    massBalanceCumulativeTopAdvOut[mbRegionI] +=
        massBalanceDt*mbLocalTopAdvOut;
    massBalanceCumulativeTopCorrectedOut[mbRegionI] +=
        massBalanceDt*mbLocalTopCorrectedOut;
    massBalanceCumulativeMacroBirthMass[mbRegionI] +=
        massBalanceDt*mbLocalMacroBirthMass;
    massBalanceCumulativeMacroBirthVolume[mbRegionI] +=
        massBalanceDt*mbLocalMacroBirthVolume;

    const scalar mbDeltaMass =
        mbLocalMass - massBalanceLastMass[mbRegionI];
    const scalar mbNetFieldMass =
        mbLocalMass - massBalanceInitialMass[mbRegionI];
    massBalanceLastMass[mbRegionI] = mbLocalMass;

    // All reductions and accumulations occur at EVERY time step;
    // only printing is decimated.
    if
    (
        Pstream::master()
     && (
            massBalanceStepCount % massBalanceReportInterval == 0
         || runTime.writeTime()
        )
    )
    {
        Info<< "MASSCHK_STEP"
            << " region=" << mbMesh.name()
            << " step=" << massBalanceStepCount
            << " time=" << runTime.value()
            << " dt=" << massBalanceDt
            << " massKg=" << mbLocalMass
            << " deltaMassStepKg=" << mbDeltaMass
            << " deltaMassTotalKg=" << mbNetFieldMass
            << " vaporVolumeM3=" << mbLocalVaporVolume
            << " vaporMassKg=" << mbLocalVaporMass
            << " advBoundaryOutKgS=" << mbLocalAdvOut
            << " correctedBoundaryOutKgS=" << mbLocalCorrectedOut
            << " poolTopAdvOutKgS=" << mbLocalTopAdvOut
            << " poolTopCorrectedOutKgS=" << mbLocalTopCorrectedOut
            << " combinedPhaseChangeMassKgS=" << mbLocalPhaseSource
            << " resolvedPsi0KgS=" << mbLocalResolvedPsi0
            << " macroBirthVaporMassKgS=" << mbLocalMacroBirthMass
            << " macroBirthVaporVolumeM3S=" << mbLocalMacroBirthVolume
            << " cumulativeAdvOutKg="
            << massBalanceCumulativeAdvOut[mbRegionI]
            << " cumulativeCorrectedOutKg="
            << massBalanceCumulativeCorrectedOut[mbRegionI]
            << " cumulativeTopAdvOutKg="
            << massBalanceCumulativeTopAdvOut[mbRegionI]
            << " cumulativeTopCorrectedOutKg="
            << massBalanceCumulativeTopCorrectedOut[mbRegionI]
            << " cumulativeMacroBirthVaporMassKg="
            << massBalanceCumulativeMacroBirthMass[mbRegionI]
            << " cumulativeMacroBirthVaporVolumeM3="
            << massBalanceCumulativeMacroBirthVolume[mbRegionI]
            << " fieldPlusAdvOutKg="
            << mbNetFieldMass + massBalanceCumulativeAdvOut[mbRegionI]
            << " fieldPlusCorrectedOutKg="
            << mbNetFieldMass + massBalanceCumulativeCorrectedOut[mbRegionI]
            << endl;
    }
}
