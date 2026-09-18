/*---------------------------------------------------------------------------*\
    Coarse persistent-population boiling model with energy-limited detached-bubble release
\*---------------------------------------------------------------------------*/

#include "populationDetachedBubbleBirth.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvPatch.H"
#include "fvScalarMatrix.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "surfaceFields.H"

#include <cmath>

namespace Foam
{
    defineTypeNameAndDebug(populationDetachedBubbleBirth, 0);
    addToRunTimeSelectionTable
    (
        macroModel,
        populationDetachedBubbleBirth,
        components
    );
}

Foam::populationDetachedBubbleBirth::populationDetachedBubbleBirth
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& p,
    singleComponentSatProp& satModel,
    const compressibleInterPhaseTransportModel& turbModel,
    const dictionary& dict
)
:
    macroModel
    (
        typeName,
        phase1,
        phase2,
        p,
        satModel,
        turbModel,
        dict
    ),
    patch_(modelDict().get<word>("patch")),
    siteCentres_(modelDict().lookup("siteCentres")),
    energyOwnershipMode_
    (
        modelDict().lookupOrDefault<word>("energyOwnershipMode", "captureRadius")
    ),
    energyOwnershipCentres_(),
    siteOwnershipIds_(),
    TSatValue_(modelDict().lookupOrDefault<scalar>("Tsat", -1)),
    siteActivationSuperheat_
    (
        modelDict().lookup("siteActivationSuperheat")
    ),
    activationPersistence_
    (
        modelDict().lookupOrDefault<scalar>("activationPersistence", 5.0e-4)
    ),
    activationStagger_
    (
        modelDict().lookupOrDefault<scalar>("activationStagger", 5.0e-4)
    ),
    captureRadius_
    (
        modelDict().lookupOrDefault<scalar>("captureRadius", 8.0e-4)
    ),
    departureRadius_
    (
        modelDict().lookupOrDefault<scalar>("departureRadius", 5.0e-4)
    ),
    detachmentGap_
    (
        modelDict().lookupOrDefault<scalar>("detachmentGap", 1.0e-4)
    ),
    interfaceWidthCells_
    (
        modelDict().lookupOrDefault<scalar>("interfaceWidthCells", 1.0)
    ),
    targetVaporFraction_
    (
        modelDict().lookupOrDefault<scalar>("targetVaporFraction", 0.9995)
    ),
    maxAlphaVaporPerStep_
    (
        modelDict().lookupOrDefault<scalar>("maxAlphaVaporPerStep", 0.03)
    ),
    creationSteps_(modelDict().lookupOrDefault<label>("creationSteps", 30)),
    maximumCreationSteps_
    (
        modelDict().lookupOrDefault<label>("maximumCreationSteps", 300)
    ),
    wallFluxUseFraction_
    (
        modelDict().lookupOrDefault<scalar>("wallFluxUseFraction", 1.0)
    ),
    sensibleReserveFraction_
    (
        modelDict().lookupOrDefault<scalar>("sensibleReserveFraction", 0.25)
    ),
    minimumRearmDelay_
    (
        modelDict().lookupOrDefault<scalar>("minimumRearmDelay", 2.0e-3)
    ),
    shieldingRadius_
    (
        modelDict().lookupOrDefault<scalar>("shieldingRadius", 1.4e-3)
    ),
    shieldingHeight_
    (
        modelDict().lookupOrDefault<scalar>("shieldingHeight", 6.0e-3)
    ),
    shieldingAlphaVaporThreshold_
    (
        modelDict().lookupOrDefault<scalar>
        (
            "shieldingAlphaVaporThreshold",
            0.02
        )
    ),
    continuitySourceMode_
    (
        modelDict().lookupOrDefault<word>("continuitySourceMode", "vaporMass")
    ),
    writeDiagnostics_
    (
        modelDict().lookupOrDefault<Switch>("writeDiagnostics", true)
    ),
    diagnosticPrefix_
    (
        modelDict().lookupOrDefault<word>
        (
            "diagnosticPrefix",
            "fixedSiteRelease"
        )
    ),
    diagnosticsInterval_
    (
        modelDict().lookupOrDefault<label>("diagnosticsInterval", 50)
    ),
    lastUpdateTimeIndex_(-1),
    patchID_(-1),
    siteFaceLabels_(),
    sphereCells_(),
    sphereTargetVapor_(),
    shieldingCells_(),
    captureArea_(),
    targetVaporVolume_(),
    targetVaporMass_(),
    targetLatentEnergy_(),
    siteState_(),
    activationDuration_(),
    storedEnergy_(),
    releaseTime_(),
    nextActivationTime_(),
    creationStep_(),
    cycleId_(),
    createdVaporVolume_(),
    createdVaporMass_(),
    consumedLatentEnergy_(),
    activeReleaseMask_
    (
        IOobject
        (
            word(diagnosticPrefix_ + "ActiveReleaseMask"),
            phase1.mesh().time().timeName(),
            phase1.mesh(),
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        phase1.mesh(),
        dimensionedScalar("zero", dimless, 0.0),
        "zeroGradient"
    ),
    siteStateField_
    (
        IOobject
        (
            word(diagnosticPrefix_ + "SiteState"),
            phase1.mesh().time().timeName(),
            phase1.mesh(),
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        phase1.mesh(),
        dimensionedScalar("zero", dimless, 0.0),
        "zeroGradient"
    ),
    alphaBirthSource_
    (
        IOobject
        (
            word(diagnosticPrefix_ + "AlphaSource"),
            phase1.mesh().time().timeName(),
            phase1.mesh(),
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        phase1.mesh(),
        dimensionedScalar("zero", dimless/dimTime, 0.0),
        "zeroGradient"
    ),
    massBirthSource_
    (
        IOobject
        (
            word(diagnosticPrefix_ + "MassSource"),
            phase1.mesh().time().timeName(),
            phase1.mesh(),
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        phase1.mesh(),
        dimensionedScalar("zero", dimDensity/dimTime, 0.0),
        "zeroGradient"
    ),
    latentSink_
    (
        IOobject
        (
            word(diagnosticPrefix_ + "LatentSink"),
            phase1.mesh().time().timeName(),
            phase1.mesh(),
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        phase1.mesh(),
        dimensionedScalar("zero", dimPower/dimVol, 0.0),
        "zeroGradient"
    )
{
    if (energyOwnershipMode_ == "persistentVoronoi")
    {
        energyOwnershipCentres_ =
            List<vector>(modelDict().lookup("energyOwnershipCentres"));
        siteOwnershipIds_ =
            labelList(modelDict().lookup("siteOwnershipIds"));
    }

    validateControls();
    buildSiteStencils();

    const label nSites = siteCentres_.size();
    siteState_.setSize(nSites, WAITING_FOR_ACTIVATION);
    activationDuration_.setSize(nSites, 0.0);
    storedEnergy_.setSize(nSites, 0.0);
    releaseTime_.setSize(nSites, -GREAT);
    nextActivationTime_.setSize(nSites, 0.0);
    creationStep_.setSize(nSites, 0);
    cycleId_.setSize(nSites, 0);
    createdVaporVolume_.setSize(nSites, 0.0);
    createdVaporMass_.setSize(nSites, 0.0);
    consumedLatentEnergy_.setSize(nSites, 0.0);

    forAll(nextActivationTime_, siteI)
    {
        nextActivationTime_[siteI] = activationStagger_*siteI;
    }

    if (Pstream::master())
    {
        Info<< "populationDetachedBubbleBirth configured: patch=" << patch_
            << ", sites=" << nSites
            << ", departureRadius=" << departureRadius_
            << ", detachmentGap=" << detachmentGap_
            << ", captureRadius=" << captureRadius_
            << ", energyOwnershipMode=" << energyOwnershipMode_
            << ", ownershipCandidates=" << energyOwnershipCentres_.size()
            << ", activationThresholds=" << siteActivationSuperheat_.size()
            << ", targetVaporFraction=" << targetVaporFraction_
            << ", continuitySourceMode=" << continuitySourceMode_
            << ", parallel=" << Pstream::parRun()
            << endl;
    }
}

void Foam::populationDetachedBubbleBirth::validateControls() const
{
    if (siteCentres_.empty())
    {
        FatalErrorInFunction
            << "siteCentres must contain at least one physical site"
            << exit(FatalError);
    }
    if
    (
        energyOwnershipMode_ != "captureRadius"
     && energyOwnershipMode_ != "persistentVoronoi"
    )
    {
        FatalErrorInFunction
            << "energyOwnershipMode must be captureRadius or persistentVoronoi; found "
            << energyOwnershipMode_
            << exit(FatalError);
    }
    if (energyOwnershipMode_ == "persistentVoronoi")
    {
        if (energyOwnershipCentres_.empty())
        {
            FatalErrorInFunction
                << "persistentVoronoi requires non-empty energyOwnershipCentres"
                << exit(FatalError);
        }
        if (siteOwnershipIds_.size() != siteCentres_.size())
        {
            FatalErrorInFunction
                << "siteOwnershipIds must contain exactly one persistent owner "
                << "index per siteCentre; sites=" << siteCentres_.size()
                << ", ownershipIds=" << siteOwnershipIds_.size()
                << exit(FatalError);
        }
        forAll(siteOwnershipIds_, siteI)
        {
            const label ownerI = siteOwnershipIds_[siteI];
            if (ownerI < 0 || ownerI >= energyOwnershipCentres_.size())
            {
                FatalErrorInFunction
                    << "siteOwnershipIds[" << siteI << "]=" << ownerI
                    << " is outside energyOwnershipCentres range [0,"
                    << energyOwnershipCentres_.size() - 1 << "]"
                    << exit(FatalError);
            }
            if
            (
                mag(siteCentres_[siteI] - energyOwnershipCentres_[ownerI])
              > 1.0e-9
            )
            {
                FatalErrorInFunction
                    << "siteCentre " << siteI << "=" << siteCentres_[siteI]
                    << " does not match mapped persistent ownership centre "
                    << ownerI << "=" << energyOwnershipCentres_[ownerI]
                    << exit(FatalError);
            }
            for (label otherI = siteI + 1; otherI < siteOwnershipIds_.size(); ++otherI)
            {
                if (siteOwnershipIds_[otherI] == ownerI)
                {
                    FatalErrorInFunction
                        << "persistent ownership index " << ownerI
                        << " is mapped to more than one active site"
                        << exit(FatalError);
                }
            }
        }
    }
    if (siteActivationSuperheat_.size() != siteCentres_.size())
    {
        FatalErrorInFunction
            << "siteActivationSuperheat must contain exactly one threshold per "
            << "siteCentre; sites=" << siteCentres_.size()
            << ", thresholds=" << siteActivationSuperheat_.size()
            << exit(FatalError);
    }
    forAll(siteActivationSuperheat_, siteI)
    {
        if (siteActivationSuperheat_[siteI] < 0)
        {
            FatalErrorInFunction
                << "siteActivationSuperheat[" << siteI
                << "] must be non-negative; found "
                << siteActivationSuperheat_[siteI]
                << exit(FatalError);
        }
    }
    if (activationPersistence_ < 0)
    {
        FatalErrorInFunction
            << "activationPersistence must be non-negative"
            << exit(FatalError);
    }
    if (activationStagger_ < 0)
    {
        FatalErrorInFunction
            << "activationStagger must be non-negative"
            << exit(FatalError);
    }
    if
    (
        captureRadius_ <= SMALL
     || departureRadius_ <= SMALL
     || detachmentGap_ < 0
     || interfaceWidthCells_ <= SMALL
    )
    {
        FatalErrorInFunction
            << "captureRadius/departureRadius/interfaceWidthCells must be positive "
            << "and detachmentGap must be non-negative"
            << exit(FatalError);
    }
    if (targetVaporFraction_ <= 0 || targetVaporFraction_ >= 1)
    {
        FatalErrorInFunction
            << "targetVaporFraction must be in (0,1)"
            << exit(FatalError);
    }
    if (maxAlphaVaporPerStep_ <= 0 || maxAlphaVaporPerStep_ > 1)
    {
        FatalErrorInFunction
            << "maxAlphaVaporPerStep must be in (0,1]"
            << exit(FatalError);
    }
    if (creationSteps_ < 1 || maximumCreationSteps_ < creationSteps_)
    {
        FatalErrorInFunction
            << "maximumCreationSteps must be >= creationSteps >= 1"
            << exit(FatalError);
    }
    if
    (
        wallFluxUseFraction_ < 0 || wallFluxUseFraction_ > 1
     || sensibleReserveFraction_ < 0 || sensibleReserveFraction_ > 1
    )
    {
        FatalErrorInFunction
            << "wallFluxUseFraction and sensibleReserveFraction must be in [0,1]"
            << exit(FatalError);
    }
    if
    (
        minimumRearmDelay_ < 0
     || shieldingRadius_ <= SMALL
     || shieldingHeight_ <= SMALL
     || shieldingAlphaVaporThreshold_ < 0
     || shieldingAlphaVaporThreshold_ > 1
    )
    {
        FatalErrorInFunction
            << "Invalid rearm/shielding controls"
            << exit(FatalError);
    }
    if
    (
        continuitySourceMode_ != "vaporMass"
     && continuitySourceMode_ != "none"
    )
    {
        FatalErrorInFunction
            << "continuitySourceMode must be vaporMass or none"
            << exit(FatalError);
    }
    if (diagnosticsInterval_ < 1)
    {
        FatalErrorInFunction
            << "diagnosticsInterval must be >= 1"
            << exit(FatalError);
    }
}

Foam::scalar Foam::populationDetachedBubbleBirth::TSat(const label celli) const
{
    return TSatValue_ > 0 ? TSatValue_ : satModel_.TSat()[celli];
}

const char* Foam::populationDetachedBubbleBirth::stateName
(
    const label state
) const
{
    if (state == WAITING_FOR_ACTIVATION) return "WAITING_FOR_ACTIVATION";
    if (state == CHARGING_DEPARTURE_BUBBLE) return "CHARGING_DEPARTURE_BUBBLE";
    if (state == CREATING_DETACHED_BUBBLE) return "CREATING_DETACHED_BUBBLE";
    if (state == WAITING_FOR_CLEARANCE) return "WAITING_FOR_CLEARANCE";
    return "UNKNOWN";
}

void Foam::populationDetachedBubbleBirth::buildSiteStencils()
{
    const fvMesh& mesh = phase1_.mesh();
    patchID_ = mesh.boundaryMesh().findPatchID(patch_);
    if (patchID_ < 0)
    {
        FatalErrorInFunction
            << "Cannot find configured heater patch " << patch_
            << " in region " << mesh.name()
            << exit(FatalError);
    }

    const fvPatch& patch = mesh.boundary()[patchID_];
    const vectorField& faceCentres = patch.Cf();
    const vectorField& faceAreaVectors = patch.Sf();

    labelList faceOwnership(faceCentres.size(), -1);
    if (energyOwnershipMode_ == "persistentVoronoi")
    {
        forAll(faceCentres, faceI)
        {
            scalar nearestDistanceSqr = GREAT;
            label nearestOwner = -1;
            forAll(energyOwnershipCentres_, ownerI)
            {
                const vector d =
                    faceCentres[faceI] - energyOwnershipCentres_[ownerI];
                const scalar radialDistanceSqr = sqr(d.x()) + sqr(d.y());
                if (radialDistanceSqr < nearestDistanceSqr)
                {
                    nearestDistanceSqr = radialDistanceSqr;
                    nearestOwner = ownerI;
                }
            }
            faceOwnership[faceI] = nearestOwner;
        }
    }

    tmp<volScalarField> tRhoV = phase2_.thermo().rho();
    const volScalarField& rhoV = tRhoV();

    const label nSites = siteCentres_.size();
    siteFaceLabels_.setSize(nSites);
    sphereCells_.setSize(nSites);
    sphereTargetVapor_.setSize(nSites);
    shieldingCells_.setSize(nSites);
    captureArea_.setSize(nSites, 0.0);
    targetVaporVolume_.setSize(nSites, 0.0);
    targetVaporMass_.setSize(nSites, 0.0);
    targetLatentEnergy_.setSize(nSites, 0.0);

    forAll(siteCentres_, siteI)
    {
        const vector& site = siteCentres_[siteI];
        DynamicList<label> localFaces;
        scalar localCaptureArea = 0;
        forAll(faceCentres, faceI)
        {
            bool ownsFace = false;
            if (energyOwnershipMode_ == "persistentVoronoi")
            {
                ownsFace =
                    faceOwnership[faceI] == siteOwnershipIds_[siteI];
            }
            else
            {
                const vector d = faceCentres[faceI] - site;
                const scalar radial =
                    Foam::sqrt(sqr(d.x()) + sqr(d.y()));
                ownsFace = radial <= captureRadius_ + SMALL;
            }

            if (ownsFace)
            {
                localFaces.append(faceI);
                localCaptureArea += mag(faceAreaVectors[faceI]);
            }
        }
        siteFaceLabels_[siteI].transfer(localFaces);
        reduce(localCaptureArea, sumOp<scalar>());
        captureArea_[siteI] = localCaptureArea;

        if
        (
            energyOwnershipMode_ == "persistentVoronoi"
         && Pstream::master()
        )
        {
            Info<< "POPULATION_ENERGY_OWNERSHIP"
                << " site=" << siteI
                << " persistentSite=" << siteOwnershipIds_[siteI]
                << " centre=" << site
                << " area=" << captureArea_[siteI]
                << endl;
        }

        const vector bubbleCentre =
            site + vector(0, 0, departureRadius_ + detachmentGap_);
        DynamicList<label> localSphereCells;
        DynamicList<scalar> localTargets;
        scalar localVolume = 0;
        scalar localMass = 0;
        scalar localLatent = 0;
        label localCoreCells = 0;

        forAll(mesh.C(), celli)
        {
            const scalar cellSize = std::cbrt(max(mesh.V()[celli], VSMALL));
            const scalar interfaceWidth =
                max(interfaceWidthCells_*cellSize, scalar(SMALL));
            const scalar signedDistance =
                departureRadius_ - mag(mesh.C()[celli] - bubbleCentre);
            const scalar profile = min
            (
                max
                (
                    scalar(0.5) + signedDistance/interfaceWidth,
                    scalar(0)
                ),
                scalar(1)
            );
            const scalar target = min(targetVaporFraction_, profile);
            if (target > SMALL)
            {
                localSphereCells.append(celli);
                localTargets.append(target);
                const scalar vaporVolume = target*mesh.V()[celli];
                const scalar vaporMass = rhoV[celli]*vaporVolume;
                localVolume += vaporVolume;
                localMass += vaporMass;
                localLatent += vaporMass*satModel_.L()[celli];
                localCoreCells += target > 0.95;
            }
        }
        sphereCells_[siteI].transfer(localSphereCells);
        sphereTargetVapor_[siteI].setSize(localTargets.size());
        forAll(sphereTargetVapor_[siteI], targetI)
        {
            sphereTargetVapor_[siteI][targetI] = localTargets[targetI];
        }

        reduce(localVolume, sumOp<scalar>());
        reduce(localMass, sumOp<scalar>());
        reduce(localLatent, sumOp<scalar>());
        reduce(localCoreCells, sumOp<label>());
        targetVaporVolume_[siteI] = localVolume;
        targetVaporMass_[siteI] = localMass;
        targetLatentEnergy_[siteI] = localLatent;

        DynamicList<label> localShieldCells;
        forAll(mesh.C(), celli)
        {
            const vector d = mesh.C()[celli] - site;
            const scalar radial = Foam::sqrt(sqr(d.x()) + sqr(d.y()));
            if
            (
                d.z() >= -SMALL
             && d.z() <= shieldingHeight_ + SMALL
             && radial <= shieldingRadius_ + SMALL
            )
            {
                localShieldCells.append(celli);
            }
        }
        shieldingCells_[siteI].transfer(localShieldCells);

        label globalFaceCount = siteFaceLabels_[siteI].size();
        label globalSphereCount = sphereCells_[siteI].size();
        label globalShieldCount = shieldingCells_[siteI].size();
        reduce(globalFaceCount, sumOp<label>());
        reduce(globalSphereCount, sumOp<label>());
        reduce(globalShieldCount, sumOp<label>());

        if
        (
            globalFaceCount == 0
         || globalSphereCount == 0
         || globalShieldCount == 0
         || localCaptureArea <= SMALL
         || localVolume <= SMALL
         || localLatent <= SMALL
         || localCoreCells == 0
        )
        {
            FatalErrorInFunction
                << "Site " << siteI << " at " << site
                << " has an invalid capture, release, or shielding stencil"
                << exit(FatalError);
        }
    }
}

void Foam::populationDetachedBubbleBirth::resetSite
(
    const label siteI,
    const scalar timeValue
)
{
    siteState_[siteI] = WAITING_FOR_ACTIVATION;
    activationDuration_[siteI] = 0;
    storedEnergy_[siteI] = 0;
    creationStep_[siteI] = 0;
    createdVaporVolume_[siteI] = 0;
    createdVaporMass_[siteI] = 0;
    consumedLatentEnergy_[siteI] = 0;
    releaseTime_[siteI] = -GREAT;
    nextActivationTime_[siteI] = timeValue + minimumRearmDelay_;
}

void Foam::populationDetachedBubbleBirth::updateSources()
{
    const fvMesh& mesh = phase1_.mesh();
    const label timeIndex = mesh.time().timeIndex();
    if (timeIndex == lastUpdateTimeIndex_)
    {
        return;
    }
    lastUpdateTimeIndex_ = timeIndex;

    activeReleaseMask_ *= scalar(0);
    siteStateField_ *= scalar(0);
    alphaBirthSource_ *= scalar(0);
    massBirthSource_ *= scalar(0);
    latentSink_ *= scalar(0);

    const scalar timeValue = mesh.time().value();
    const scalar deltaT = max(mesh.time().deltaTValue(), scalar(VSMALL));
    const volScalarField& liquidTemperature = phase1_.thermo().T();
    tmp<volScalarField> tRhoL = phase1_.thermo().rho();
    const volScalarField& rhoL = tRhoL();
    tmp<volScalarField> tCpL = phase1_.thermo().Cp();
    const volScalarField& cpL = tCpL();
    tmp<volScalarField> tRhoV = phase2_.thermo().rho();
    const volScalarField& rhoV = tRhoV();

    const fvPatch& patch = mesh.boundary()[patchID_];
    const labelUList& faceCells = patch.faceCells();
    const vectorField& faceAreaVectors = patch.Sf();
    const scalarField& wallTemperature =
        liquidTemperature.boundaryField()[patchID_];
    tmp<volScalarField> tAlphat = turbModel_.alphat();
    const scalarField& alphatPatch = tAlphat().boundaryField()[patchID_];
    tmp<scalarField> tKappaPatch =
        phase1_.thermo().kappaEff(alphatPatch, patchID_);
    const scalarField& kappaPatch = tKappaPatch();
    tmp<scalarField> tSnGrad =
        liquidTemperature.boundaryField()[patchID_].snGrad();
    const scalarField& snGrad = tSnGrad();

    forAll(siteCentres_, siteI)
    {
        scalar localSuperheatArea = 0;
        scalar localArea = 0;
        scalar localWallEnergy = 0;
        const labelList& faces = siteFaceLabels_[siteI];
        forAll(faces, faceListI)
        {
            const label faceI = faces[faceListI];
            const label celli = faceCells[faceI];
            const scalar area = mag(faceAreaVectors[faceI]);
            const scalar superheat =
                max(wallTemperature[faceI] - TSat(celli), scalar(0));
            localSuperheatArea += superheat*area;
            localArea += area;
            localWallEnergy +=
                wallFluxUseFraction_
               *max(kappaPatch[faceI]*snGrad[faceI]*area, scalar(0))
               *deltaT;
        }
        reduce(localSuperheatArea, sumOp<scalar>());
        reduce(localArea, sumOp<scalar>());
        reduce(localWallEnergy, sumOp<scalar>());
        const scalar meanWallSuperheat =
            localArea > SMALL ? localSuperheatArea/localArea : scalar(0);

        scalar localShieldVapor = 0;
        const labelList& shieldCells = shieldingCells_[siteI];
        forAll(shieldCells, shieldI)
        {
            const label celli = shieldCells[shieldI];
            localShieldVapor = max
            (
                localShieldVapor,
                scalar(1) - max(min(phase1_[celli], scalar(1)), scalar(0))
            );
        }
        reduce(localShieldVapor, maxOp<scalar>());
        const bool shielded =
            localShieldVapor + SMALL >= shieldingAlphaVaporThreshold_;

        scalar localSensibleEnergy = 0;
        const labelList& releaseCells = sphereCells_[siteI];
        forAll(releaseCells, releaseI)
        {
            const label celli = releaseCells[releaseI];
            localSensibleEnergy +=
                sensibleReserveFraction_
               *phase1_[celli]
               *rhoL[celli]
               *cpL[celli]
               *mesh.V()[celli]
               *max(liquidTemperature[celli] - TSat(celli), scalar(0));
        }
        reduce(localSensibleEnergy, sumOp<scalar>());

        const label oldState = siteState_[siteI];

        if (siteState_[siteI] == WAITING_FOR_ACTIVATION)
        {
            const bool activationCandidate =
                timeValue + SMALL >= nextActivationTime_[siteI]
             && !shielded
             && meanWallSuperheat + SMALL
                >= siteActivationSuperheat_[siteI];
            if (activationCandidate)
            {
                activationDuration_[siteI] += deltaT;
            }
            else
            {
                activationDuration_[siteI] = 0;
            }

            if
            (
                activationDuration_[siteI] + SMALL
             >= activationPersistence_
            )
            {
                siteState_[siteI] = CHARGING_DEPARTURE_BUBBLE;
                storedEnergy_[siteI] = localSensibleEnergy;
                activationDuration_[siteI] = 0;
                cycleId_[siteI]++;
            }
        }
        else if (siteState_[siteI] == CHARGING_DEPARTURE_BUBBLE)
        {
            if (shielded)
            {
                resetSite(siteI, timeValue);
            }
            else
            {
                storedEnergy_[siteI] += localWallEnergy;
                if
                (
                    storedEnergy_[siteI] + VSMALL
                 >= targetLatentEnergy_[siteI]
                )
                {
                    siteState_[siteI] = CREATING_DETACHED_BUBBLE;
                    creationStep_[siteI] = 0;
                }
            }
        }
        else if (siteState_[siteI] == CREATING_DETACHED_BUBBLE)
        {
            storedEnergy_[siteI] += localWallEnergy;
            creationStep_[siteI]++;

            scalar localRemainingVolume = 0;
            scalar localRequestedEnergy = 0;
            scalarField requestedAlpha(releaseCells.size(), 0.0);
            const label remainingSteps = max
            (
                creationSteps_ - creationStep_[siteI] + 1,
                label(1)
            );

            forAll(releaseCells, releaseI)
            {
                const label celli = releaseCells[releaseI];
                const scalar target = sphereTargetVapor_[siteI][releaseI];
                const scalar alphaLiquid =
                    max(min(phase1_[celli], scalar(1)), scalar(0));
                const scalar alphaVapor = scalar(1) - alphaLiquid;
                const scalar remaining = max(target - alphaVapor, scalar(0));
                localRemainingVolume += remaining*mesh.V()[celli];
                if (remaining <= SMALL || alphaLiquid <= SMALL)
                {
                    continue;
                }
                const scalar alphaStep = min
                (
                    remaining/scalar(remainingSteps),
                    min(maxAlphaVaporPerStep_, alphaLiquid)
                );
                requestedAlpha[releaseI] = alphaStep;
                localRequestedEnergy +=
                    rhoV[celli]
                   *alphaStep
                   *mesh.V()[celli]
                   *satModel_.L()[celli];
            }

            reduce(localRemainingVolume, sumOp<scalar>());
            reduce(localRequestedEnergy, sumOp<scalar>());

            const scalar volumeTolerance = max
            (
                SMALL,
                1.0e-6*targetVaporVolume_[siteI]
            );
            const scalar energyTolerance = max
            (
                SMALL,
                1.0e-6*targetLatentEnergy_[siteI]
            );
            const scalar remainingLatentBudget = max
            (
                targetLatentEnergy_[siteI]
              - consumedLatentEnergy_[siteI],
                scalar(0)
            );
            const bool stencilComplete =
                localRemainingVolume <= volumeTolerance;
            const bool latentBudgetComplete =
                remainingLatentBudget <= energyTolerance;

            if (stencilComplete || latentBudgetComplete)
            {
                siteState_[siteI] = WAITING_FOR_CLEARANCE;
                releaseTime_[siteI] = timeValue;
                if (Pstream::master())
                {
                    Info<< "POPULATION_SITE_RELEASE"
                        << " site=" << siteI
                        << " cycle=" << cycleId_[siteI]
                        << " time=" << timeValue
                        << " centre=" << siteCentres_[siteI]
                        << " departureRadius=" << departureRadius_
                        << " completion="
                        << (latentBudgetComplete ? "latentBudget" : "stencil")
                        << " createdVaporVolume="
                        << createdVaporVolume_[siteI]
                        << " createdVaporMass=" << createdVaporMass_[siteI]
                        << " consumedLatentEnergy="
                        << consumedLatentEnergy_[siteI]
                        << " targetLatentEnergy="
                        << targetLatentEnergy_[siteI]
                        << endl;
                }
            }
            else
            {
                const scalar availableEnergyScale =
                    localRequestedEnergy > SMALL
                  ? min
                    (
                        scalar(1),
                        storedEnergy_[siteI]/localRequestedEnergy
                    )
                  : scalar(0);
                const scalar budgetEnergyScale =
                    localRequestedEnergy > SMALL
                  ? min
                    (
                        scalar(1),
                        remainingLatentBudget/localRequestedEnergy
                    )
                  : scalar(0);
                const scalar energyScale = min
                (
                    availableEnergyScale,
                    budgetEnergyScale
                );

                scalar localCreatedVolume = 0;
                scalar localCreatedMass = 0;
                scalar localConsumedEnergy = 0;
                forAll(releaseCells, releaseI)
                {
                    const scalar alphaStep =
                        energyScale*requestedAlpha[releaseI];
                    if (alphaStep <= SMALL)
                    {
                        continue;
                    }
                    const label celli = releaseCells[releaseI];
                    const scalar alphaRate = alphaStep/deltaT;
                    const scalar vaporMassRate = rhoV[celli]*alphaRate;
                    const scalar latentRate =
                        vaporMassRate*satModel_.L()[celli];

                    alphaBirthSource_[celli] += alphaRate;
                    if (continuitySourceMode_ == "vaporMass")
                    {
                        massBirthSource_[celli] += vaporMassRate;
                    }
                    latentSink_[celli] += latentRate;
                    activeReleaseMask_[celli] = scalar(1);

                    const scalar createdVolume =
                        alphaStep*mesh.V()[celli];
                    const scalar createdMass =
                        rhoV[celli]*createdVolume;
                    const scalar consumedEnergy =
                        createdMass*satModel_.L()[celli];
                    localCreatedVolume += createdVolume;
                    localCreatedMass += createdMass;
                    localConsumedEnergy += consumedEnergy;
                }

                reduce(localCreatedVolume, sumOp<scalar>());
                reduce(localCreatedMass, sumOp<scalar>());
                reduce(localConsumedEnergy, sumOp<scalar>());
                createdVaporVolume_[siteI] += localCreatedVolume;
                createdVaporMass_[siteI] += localCreatedMass;
                consumedLatentEnergy_[siteI] += localConsumedEnergy;
                storedEnergy_[siteI] = max
                (
                    storedEnergy_[siteI] - localConsumedEnergy,
                    scalar(0)
                );

                const scalar postSourceRemainingLatentBudget = max
                (
                    targetLatentEnergy_[siteI]
                  - consumedLatentEnergy_[siteI],
                    scalar(0)
                );
                if
                (
                    creationStep_[siteI] >= maximumCreationSteps_
                 && postSourceRemainingLatentBudget > energyTolerance
                )
                {
                    FatalErrorInFunction
                        << "Site " << siteI
                        << " failed to create its detached bubble after "
                        << creationStep_[siteI] << " source updates"
                        << "; remainingVolume=" << localRemainingVolume
                        << "; createdVolume=" << createdVaporVolume_[siteI]
                        << "; consumedLatentEnergy="
                        << consumedLatentEnergy_[siteI]
                        << "; targetLatentEnergy="
                        << targetLatentEnergy_[siteI]
                        << "; remainingLatentBudget="
                        << postSourceRemainingLatentBudget
                        << exit(FatalError);
                }
            }
        }
        else if (siteState_[siteI] == WAITING_FOR_CLEARANCE)
        {
            if
            (
                !shielded
             && timeValue - releaseTime_[siteI] + SMALL
                >= minimumRearmDelay_
            )
            {
                resetSite(siteI, timeValue);
            }
        }

        const scalar stateValue = siteState_[siteI];
        forAll(releaseCells, releaseI)
        {
            const label celli = releaseCells[releaseI];
            siteStateField_[celli] = max
            (
                siteStateField_[celli],
                stateValue
            );
        }

        if
        (
            writeDiagnostics_
         && Pstream::master()
         &&
            (
                oldState != siteState_[siteI]
             || timeIndex % diagnosticsInterval_ == 0
             || mesh.time().writeTime()
            )
        )
        {
            Info<< "POPULATION_SITE_STATE"
                << " site=" << siteI
                << " cycle=" << cycleId_[siteI]
                << " time=" << timeValue
                << " state=" << stateName(siteState_[siteI])
                << " meanWallSuperheat=" << meanWallSuperheat
                << " activationSuperheat="
                << siteActivationSuperheat_[siteI]
                << " shieldVapor=" << localShieldVapor
                << " captureArea=" << captureArea_[siteI]
                << " storedEnergy=" << storedEnergy_[siteI]
                << " targetEnergy=" << targetLatentEnergy_[siteI]
                << " createdVolume=" << createdVaporVolume_[siteI]
                << " targetVolume=" << targetVaporVolume_[siteI]
                << endl;
        }
    }

    activeReleaseMask_.correctBoundaryConditions();
    siteStateField_.correctBoundaryConditions();
    alphaBirthSource_.correctBoundaryConditions();
    massBirthSource_.correctBoundaryConditions();
    latentSink_.correctBoundaryConditions();
}

void Foam::populationDetachedBubbleBirth::TSource1
(
    fvScalarMatrix& T1Eqn
)
{
    updateSources();
    T1Eqn.source() -=
        latentSink_.internalField()*phase1_.mesh().V();
}

void Foam::populationDetachedBubbleBirth::TSource2
(
    fvScalarMatrix& T2Eqn
)
{}

void Foam::populationDetachedBubbleBirth::energySource(volScalarField& Q)
{}

void Foam::populationDetachedBubbleBirth::energySource1(volScalarField& q1)
{}

void Foam::populationDetachedBubbleBirth::energySource2(volScalarField& q2)
{}

void Foam::populationDetachedBubbleBirth::massSource
(
    volScalarField& rhoSource
)
{
    updateSources();
    rhoSource += massBirthSource_;
}

void Foam::populationDetachedBubbleBirth::alphaSource
(
    volScalarField& rhoSource
)
{
    updateSources();
    rhoSource += alphaBirthSource_;
}

} // End namespace Foam

// ************************************************************************* //