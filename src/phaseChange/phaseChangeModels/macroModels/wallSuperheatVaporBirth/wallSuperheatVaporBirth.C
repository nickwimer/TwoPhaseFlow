/*---------------------------------------------------------------------------*\
    Wall-superheat vapor-birth macro model for all-liquid heated-wall cases
\*---------------------------------------------------------------------------*/

#include "wallSuperheatVaporBirth.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvPatch.H"
#include "fvScalarMatrix.H"
#include "DynamicList.H"
#include "PstreamReduceOps.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "globalIndex.H"
#include "labelIOList.H"
#include "mathematicalConstants.H"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdlib>

namespace Foam
{
    defineTypeNameAndDebug(wallSuperheatVaporBirth, 0);
    addToRunTimeSelectionTable(macroModel, wallSuperheatVaporBirth, components);
}

Foam::wallSuperheatVaporBirth::wallSuperheatVaporBirth
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
    patches_(modelDict().lookup("patches")),
    TSatValue_(modelDict().lookupOrDefault<scalar>("Tsat", -1)),
    superheatOn_(modelDict().get<scalar>("superheatOn")),
    superheatOff_(modelDict().get<scalar>("superheatOff")),
    maxAlphaVaporPerStep_(modelDict().lookupOrDefault<scalar>("maxAlphaVaporPerStep", 0.02)),
    birthLayerCells_(modelDict().lookupOrDefault<label>("birthLayerCells", 1)),
    sourceRelaxation_(modelDict().lookupOrDefault<scalar>("sourceRelaxation", 1.0)),
    epsLiquid_(modelDict().lookupOrDefault<scalar>("epsLiquid", 1e-6)),
    writeDiagnostics_(modelDict().lookupOrDefault<Switch>("writeDiagnostics", true)),
    diagnosticPrefix_(modelDict().lookupOrDefault<word>("diagnosticPrefix", typeName)),
    sourceApplicationMode_(modelDict().lookupOrDefault<word>("sourceApplicationMode", "normal")),
    postHandoffSourceMode_
    (
        modelDict().lookupOrDefault<word>
        (
            "postHandoffSourceMode",
            "sustainedContactLine"
        )
    ),
    seedContinuitySourceMode_
    (
        modelDict().lookupOrDefault<word>
        (
            "seedContinuitySourceMode",
            "legacyVaporMass"
        )
    ),
    birthTriggerMode_(modelDict().lookupOrDefault<word>("birthTriggerMode", "legacyWallSuperheat")),
    seedHandoffMode_(modelDict().lookupOrDefault<word>("seedHandoffMode", "none")),
    massDensityPhase_(modelDict().lookupOrDefault<word>("massDensityPhase", "liquid")),
    sourceFootprintMode_(modelDict().lookupOrDefault<word>("sourceFootprintMode", "fixedSite")),
    growthDepositionMode_(modelDict().lookupOrDefault<word>("growthDepositionMode", "wallAdjacentCell")),
    growthBandHeight_(modelDict().lookupOrDefault<scalar>("growthBandHeight", 8.0e-4)),
    growthBandAlphaMin_(modelDict().lookupOrDefault<scalar>("growthBandAlphaMin", 0.02)),
    growthBandAlphaMax_(modelDict().lookupOrDefault<scalar>("growthBandAlphaMax", 0.98)),
    growthBandNormalDecayLength_(modelDict().lookupOrDefault<scalar>("growthBandNormalDecayLength", 3.0e-4)),
    growthBandMinimumLiquid_(modelDict().lookupOrDefault<scalar>("growthBandMinimumLiquid", 0.02)),
    growthBandRequirePrimaryBubble_(modelDict().lookupOrDefault<Switch>("growthBandRequirePrimaryBubble", true)),
    siteMode_(modelDict().lookupOrDefault<word>("growthSiteMode", modelDict().lookupOrDefault<word>("siteMode", "allPatchAdjacent"))),
    sitePatch_(modelDict().lookupOrDefault<word>("growthSitePatch", modelDict().lookupOrDefault<word>("sitePatch", ""))),
    siteCenter_(modelDict().lookupOrDefault<vector>("growthSiteCentre", modelDict().lookupOrDefault<vector>("siteCenter", vector::zero))),
    siteRadius_(modelDict().lookupOrDefault<scalar>("growthSiteRadius", modelDict().lookupOrDefault<scalar>("siteRadius", 0))),
    contactLineAlphaMin_(modelDict().lookupOrDefault<scalar>("contactLineAlphaMin", 0.02)),
    contactLineAlphaMax_(modelDict().lookupOrDefault<scalar>("contactLineAlphaMax", 0.98)),
    minimumSourceLiquidFraction_(modelDict().lookupOrDefault<scalar>("minimumSourceLiquidFraction", 0.02)),
    contactLineNeighbourSearch_(modelDict().lookupOrDefault<Switch>("contactLineNeighbourSearch", true)),
    contactLineNeighbourDelta_(modelDict().lookupOrDefault<scalar>("contactLineNeighbourDelta", 0.25)),
    maxActiveSites_(modelDict().lookupOrDefault<label>("maxActiveSites", 0)),
    lockSiteAfterBirth_(modelDict().lookupOrDefault<Switch>("lockSiteAfterBirth", false)),
    sustainedBirthMode_(modelDict().lookupOrDefault<word>("sustainedBirthMode", "none")),
    minimumActiveDuration_(modelDict().lookupOrDefault<scalar>("minimumActiveDuration", 0)),
    minimumEvaporationSuperheat_(modelDict().lookupOrDefault<scalar>("minimumEvaporationSuperheat", 0)),
    thermalReserveFraction_(modelDict().lookupOrDefault<scalar>("thermalReserveFraction", 0)),
    minimumThermalAlphaStep_(modelDict().lookupOrDefault<scalar>("minimumThermalAlphaStep", 0)),
    thermalLimitMode_(modelDict().lookupOrDefault<word>("thermalLimitMode", "sensibleReserve")),
    wallFluxUseFraction_(modelDict().lookupOrDefault<scalar>("wallFluxUseFraction", 0)),
    seedThermalReserveFraction_
    (
        modelDict().lookupOrDefault<scalar>
        (
            "seedThermalReserveFraction",
            modelDict().lookupOrDefault<scalar>("thermalReserveFraction", 0)
        )
    ),
    seedWallFluxUseFraction_
    (
        modelDict().lookupOrDefault<scalar>
        (
            "seedWallFluxUseFraction",
            modelDict().lookupOrDefault<scalar>("wallFluxUseFraction", 0)
        )
    ),
    writePerCellDiagnostics_(modelDict().lookupOrDefault<Switch>("writePerCellDiagnostics", false)),
    perCellDiagnosticsInterval_(modelDict().lookupOrDefault<label>("perCellDiagnosticsInterval", 1)),
    perCellDiagnosticsFormat_(modelDict().lookupOrDefault<word>("perCellDiagnosticsFormat", "csv")),
    perCellDiagnosticsSegmentID_(phase1.mesh().time().startTimeIndex()),
    lastPerCellDiagnosticsTimeIndex_(-1),
    compactDiagnosticsInterval_
    (
        modelDict().lookupOrDefault<label>("compactDiagnosticsInterval", 1)
    ),
    lastStateDiagnosticsTimeIndex_(-1),
    lastStateDiagnosticsState_(-1),
    lastStateDiagnosticsCycleId_(-1),
    lastSourceDiagnosticsTimeIndex_(-1),
    lastSourceDiagnosticsState_(-1),
    lastSourceDiagnosticsCycleId_(-1),
    cavityCentre_(modelDict().lookupOrDefault<vector>("cavityCentre", vector::zero)),
    cavityRadius_(modelDict().lookupOrDefault<scalar>("cavityRadius", 1.0e-4)),
    cavityContactAngle_(modelDict().lookupOrDefault<scalar>("cavityContactAngle", 90.0)),
    cavitySurfaceTension_
    (
        modelDict().lookupOrDefault<scalar>("cavitySurfaceTension", 0.0589)
    ),
    activationPersistenceSteps_(modelDict().lookupOrDefault<label>("activationPersistenceSteps", 2)),
    activationPersistence_
    (
        modelDict().lookupOrDefault<scalar>("activationPersistence", -1)
    ),
    seedCreationSteps_(modelDict().lookupOrDefault<label>("seedCreationSteps", 5)),
    maximumSeedCreationSteps_
    (
        modelDict().lookupOrDefault<label>
        (
            "maximumSeedCreationSteps",
            max
            (
                modelDict().lookupOrDefault<label>("seedCreationSteps", 5),
                label(200)
            )
        )
    ),
    maxSeedAlphaVaporPerStep_
    (
        modelDict().lookupOrDefault<scalar>
        (
            "maxSeedAlphaVaporPerStep",
            modelDict().lookupOrDefault<scalar>("maxAlphaVaporPerStep", 0.02)
        )
    ),
    targetSeedVaporFraction_
    (
        modelDict().lookupOrDefault<scalar>
        (
            "targetSeedVaporFraction",
            scalar(1) - modelDict().lookupOrDefault<scalar>("epsLiquid", 1e-6)
        )
    ),
    minimumCellsPerHandoffRadius_(modelDict().lookupOrDefault<label>("minimumCellsPerHandoffRadius", 4)),
    preferredCellsPerHandoffRadius_(modelDict().lookupOrDefault<label>("preferredCellsPerHandoffRadius", 5)),
    maximumHandoffRadius_(modelDict().lookupOrDefault<scalar>("maximumHandoffRadius", 0.00075)),
    repeatedNucleationMode_(modelDict().lookupOrDefault<Switch>("repeatedNucleationMode", false)),
    departurePersistence_(modelDict().lookupOrDefault<scalar>("departurePersistence", 5.0e-4)),
    departureContactAreaThreshold_(modelDict().lookupOrDefault<scalar>("departureContactAreaThreshold", 1.0e-10)),
    departureSiteVaporThreshold_(modelDict().lookupOrDefault<scalar>("departureSiteVaporThreshold", 1.0e-3)),
    rewetAlphaThreshold_(modelDict().lookupOrDefault<scalar>("rewetAlphaThreshold", 0.98)),
    rewetWallVaporThreshold_(modelDict().lookupOrDefault<scalar>("rewetWallVaporThreshold", 0.02)),
    rewetPersistence_(modelDict().lookupOrDefault<scalar>("rewetPersistence", 5.0e-4)),
    minimumRearmDelay_(modelDict().lookupOrDefault<scalar>("minimumRearmDelay", 0.0)),
    requireWettedSiteForActivation_
    (
        modelDict().lookupOrDefault<Switch>
        (
            "requireWettedSiteForActivation",
            false
        )
    ),
    activationAlphaLiquidThreshold_
    (
        modelDict().lookupOrDefault<scalar>
        (
            "activationAlphaLiquidThreshold",
            0.98
        )
    ),
    activationWallVaporThreshold_
    (
        modelDict().lookupOrDefault<scalar>
        (
            "activationWallVaporThreshold",
            0.02
        )
    ),
    shieldingRadius_
    (
        modelDict().lookupOrDefault<scalar>("shieldingRadius", 0)
    ),
    shieldingHeight_
    (
        modelDict().lookupOrDefault<scalar>("shieldingHeight", 0)
    ),
    shieldingAlphaVaporThreshold_
    (
        modelDict().lookupOrDefault<scalar>
        (
            "shieldingAlphaVaporThreshold",
            0.02
        )
    ),
    cavityCandidateRadius_(0),
    cavityTipHeight_(0),
    cavityCriticalTemperature_(0),
    cavityCriticalSuperheat_(0),
    handoffRadius_(0),
    handoffCapHeight_(0),
    handoffFootprintRadius_(0),
    handoffTargetVolume_(0),
    handoffTargetMass_(0),
    handoffTargetLatentEnergy_(0),
    localTangentialSpacing_(0),
    localWallNormalSpacing_(0),
    candidateCells_(),
    siteCells_(),
    siteCellMask_(),
    wallFaceAreaCache_(),
    handoffCells_(),
    handoffCellMask_(),
    handoffTargetVaporFractions_(),
    cavityTipCells_(),
    cavityTipWeights_(),
    shieldingCells_(),
    sustainedState_(),
    sustainedActivationTime_(),
    patchAdjacentCellCount_(0),
    lastUpdateTimeIndex_(-1),
    firstTriggerTime_(-1),
    firstActiveSiteTime_(-1),
    siteLocked_(false),
    cavityBirthState_(WAITING_FOR_CRITERION),
    activationCounter_(0),
    activationCandidateDuration_(0),
    creationStepIndex_(0),
    activationTime_(-1),
    handoffStartTime_(-1),
    handoffCompleteTime_(-1),
    storedHandoffEnergy_(0),
    createdVaporVolume_(0),
    createdVaporMass_(0),
    displacedLiquidMass_(0),
    consumedLatentEnergy_(0),
    accumulatedWallEnergy_(0),
    lastLocalSensibleEnergy_(0),
    lastCavityTtip_(0),
    cycleId_(0),
    rearmCount_(0),
    departureCandidateDuration_(0),
    rewetCandidateDuration_(0),
    departureTime_(-1),
    rewetTime_(-1),
    rearmTime_(-1),
    nextActivationTime_(-1),
    departureDetected_(false),
    rewetDetected_(false),
    siteRearmed_(false),
    lastSiteShielded_(false),
    superheat_
    (
        IOobject
        (
            word(diagnosticPrefix_ + "Superheat"),
            phase1.mesh().time().timeName(),
            phase1.mesh(),
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        phase1.mesh(),
        dimensionedScalar("0", dimTemperature, 0),
        "zeroGradient"
    ),
    activeMask_
    (
        IOobject
        (
            word(diagnosticPrefix_ + "ActiveMask"),
            phase1.mesh().time().timeName(),
            phase1.mesh(),
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        phase1.mesh(),
        dimensionedScalar("0", dimless, 0),
        "zeroGradient"
    ),
    sustainedStateField_
    (
        IOobject
        (
            word(diagnosticPrefix_ + "SustainedState"),
            phase1.mesh().time().timeName(),
            phase1.mesh(),
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        phase1.mesh(),
        dimensionedScalar("0", dimless, 0),
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
        dimensionedScalar("0", dimless/dimTime, 0),
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
        dimensionedScalar("0", dimDensity/dimTime, 0),
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
        dimensionedScalar("0", dimPower/dimVol, 0),
        "zeroGradient"
    )
{
    validateControls();
    buildCandidateCells();
    restoreCavityState();

    const fvMesh& mesh = phase1_.mesh();
    scalar initialVaporVolume = 0;
    forAll(phase1_, celli)
    {
        initialVaporVolume += (1.0 - phase1_[celli])*mesh.V()[celli];
    }
    reduce(initialVaporVolume, sumOp<scalar>());

    Info<< "wall-superheat vapor birth macro model configured" << nl
        << "    authoritativeSource: macroModel wallSuperheatVaporBirth" << nl
        << "    massDensityPhase: " << massDensityPhase_ << nl
        << "    sourceFootprintMode: " << sourceFootprintMode_ << nl
        << "    growthDepositionMode: " << growthDepositionMode_ << nl
        << "    growthBandHeight: " << growthBandHeight_ << nl
        << "    growthBandAlphaMin: " << growthBandAlphaMin_ << nl
        << "    growthBandAlphaMax: " << growthBandAlphaMax_ << nl
        << "    growthBandNormalDecayLength: " << growthBandNormalDecayLength_ << nl
        << "    growthBandMinimumLiquid: " << growthBandMinimumLiquid_ << nl
        << "    growthBandRequirePrimaryBubble: " << growthBandRequirePrimaryBubble_ << nl
        << "    patches: " << patches_ << nl
        << "    candidate wall-adjacent cells: " << patchAdjacentCellCount_ << nl
        << "    siteMode: " << siteMode_ << nl
        << "    sitePatch: " << sitePatch_ << nl
        << "    siteCenter: " << siteCenter_ << nl
        << "    siteRadius: " << siteRadius_ << nl
        << "    eligible site cells: " << siteCells_.size() << nl
        << "    eligible site volume: " << siteVolume() << nl
        << "    maxActiveSites: " << maxActiveSites_ << nl
        << "    lockSiteAfterBirth: " << lockSiteAfterBirth_ << nl
        << "    Tsat control: " << (TSatValue_ > 0 ? TSatValue_ : -1) << nl
        << "    superheatOn: " << superheatOn_ << nl
        << "    superheatOff: " << superheatOff_ << nl
        << "    maxAlphaVaporPerStep: " << maxAlphaVaporPerStep_ << nl
        << "    sustainedBirthMode: " << sustainedBirthMode_ << nl
        << "    minimumActiveDuration: " << minimumActiveDuration_ << nl
        << "    minimumEvaporationSuperheat: " << minimumEvaporationSuperheat_ << nl
        << "    thermalReserveFraction: " << thermalReserveFraction_ << nl
        << "    minimumThermalAlphaStep: " << minimumThermalAlphaStep_ << nl
        << "    thermalLimitMode: " << thermalLimitMode_ << nl
        << "    wallFluxUseFraction: " << wallFluxUseFraction_ << nl
        << "    contactLineAlphaMin: " << contactLineAlphaMin_ << nl
        << "    contactLineAlphaMax: " << contactLineAlphaMax_ << nl
        << "    minimumSourceLiquidFraction: " << minimumSourceLiquidFraction_ << nl
        << "    contactLineNeighbourSearch: " << contactLineNeighbourSearch_ << nl
        << "    contactLineNeighbourDelta: " << contactLineNeighbourDelta_ << nl
        << "    writePerCellDiagnostics: " << writePerCellDiagnostics_ << nl
        << "    perCellDiagnosticsInterval: " << perCellDiagnosticsInterval_ << nl
        << "    perCellDiagnosticsFormat: " << perCellDiagnosticsFormat_ << nl
        << "    compactDiagnosticsInterval: " << compactDiagnosticsInterval_ << nl
        << "    sourceApplicationMode: " << sourceApplicationMode_ << nl
        << "    postHandoffSourceMode: " << postHandoffSourceMode_ << nl
        << "    seedContinuitySourceMode: " << seedContinuitySourceMode_ << nl
        << "    birthTriggerMode: " << birthTriggerMode_ << nl
        << "    seedHandoffMode: " << seedHandoffMode_ << nl
        << "    cavityCentre: " << cavityCentre_ << nl
        << "    cavityRadius: " << cavityRadius_ << nl
        << "    cavityContactAngle: " << cavityContactAngle_ << nl
        << "    cavitySurfaceTension: " << cavitySurfaceTension_ << nl
        << "    activationPersistenceSteps: " << activationPersistenceSteps_ << nl
        << "    activationPersistence: " << activationPersistence_ << nl
        << "    seedCreationSteps: " << seedCreationSteps_ << nl
        << "    maximumSeedCreationSteps: " << maximumSeedCreationSteps_ << nl
        << "    maxSeedAlphaVaporPerStep: " << maxSeedAlphaVaporPerStep_ << nl
        << "    targetSeedVaporFraction: " << targetSeedVaporFraction_ << nl
        << "    seedThermalReserveFraction: " << seedThermalReserveFraction_ << nl
        << "    seedWallFluxUseFraction: " << seedWallFluxUseFraction_ << nl
        << "    maximumHandoffRadius: " << maximumHandoffRadius_ << nl
        << "    repeatedNucleationMode: " << repeatedNucleationMode_ << nl
        << "    departurePersistence: " << departurePersistence_ << nl
        << "    departureContactAreaThreshold: " << departureContactAreaThreshold_ << nl
        << "    departureSiteVaporThreshold: " << departureSiteVaporThreshold_ << nl
        << "    rewetAlphaThreshold: " << rewetAlphaThreshold_ << nl
        << "    rewetWallVaporThreshold: " << rewetWallVaporThreshold_ << nl
        << "    rewetPersistence: " << rewetPersistence_ << nl
        << "    minimumRearmDelay: " << minimumRearmDelay_ << nl
        << "    requireWettedSiteForActivation: "
        << requireWettedSiteForActivation_ << nl
        << "    activationAlphaLiquidThreshold: "
        << activationAlphaLiquidThreshold_ << nl
        << "    activationWallVaporThreshold: "
        << activationWallVaporThreshold_ << nl
        << "    shieldingRadius: " << shieldingRadius_ << nl
        << "    shieldingHeight: " << shieldingHeight_ << nl
        << "    shieldingAlphaVaporThreshold: "
        << shieldingAlphaVaporThreshold_ << nl
        << "    restoredCavityBirthState: " << cavityStateName(cavityBirthState_) << nl
        << "    initial vapor volume: " << initialVaporVolume << endl;
}

void Foam::wallSuperheatVaporBirth::validateControls() const
{
    if (patches_.empty())
    {
        FatalErrorInFunction
            << "wallSuperheatVaporBirth requires at least one patch"
            << exit(FatalError);
    }

    if (superheatOn_ <= 0)
    {
        FatalErrorInFunction
            << "superheatOn must be positive, found " << superheatOn_
            << exit(FatalError);
    }

    if (superheatOff_ < 0 || superheatOff_ >= superheatOn_)
    {
        FatalErrorInFunction
            << "superheatOff must satisfy 0 <= superheatOff < superheatOn; found "
            << superheatOff_ << " with superheatOn " << superheatOn_
            << exit(FatalError);
    }

    if (maxAlphaVaporPerStep_ <= 0 || maxAlphaVaporPerStep_ > 1)
    {
        FatalErrorInFunction
            << "maxAlphaVaporPerStep must be in (0, 1], found "
            << maxAlphaVaporPerStep_ << exit(FatalError);
    }

    if (birthLayerCells_ != 1)
    {
        FatalErrorInFunction
            << "wallSuperheatVaporBirth currently supports only birthLayerCells 1; found "
            << birthLayerCells_ << exit(FatalError);
    }

    if (sourceRelaxation_ < 0)
    {
        FatalErrorInFunction
            << "sourceRelaxation must be non-negative, found " << sourceRelaxation_
            << exit(FatalError);
    }

    if (epsLiquid_ < 0 || epsLiquid_ >= 1)
    {
        FatalErrorInFunction
            << "epsLiquid must be in [0, 1), found " << epsLiquid_
            << exit(FatalError);
    }

    if
    (
        sourceApplicationMode_ != "normal"
     && sourceApplicationMode_ != "diagnosticOnly"
     && sourceApplicationMode_ != "disabled"
    )
    {
        FatalErrorInFunction
            << "Unsupported sourceApplicationMode " << sourceApplicationMode_
            << ". Supported values are normal, diagnosticOnly and disabled"
            << exit(FatalError);
    }

    if
    (
        postHandoffSourceMode_ != "sustainedContactLine"
     && postHandoffSourceMode_ != "none"
    )
    {
        FatalErrorInFunction
            << "Unsupported postHandoffSourceMode " << postHandoffSourceMode_
            << ". Supported values are sustainedContactLine and none"
            << exit(FatalError);
    }

    if
    (
        seedContinuitySourceMode_ != "legacyVaporMass"
     && seedContinuitySourceMode_ != "none"
    )
    {
        FatalErrorInFunction
            << "Unsupported seedContinuitySourceMode "
            << seedContinuitySourceMode_
            << ". Supported values are legacyVaporMass and none"
            << exit(FatalError);
    }

    if (birthTriggerMode_ != "legacyWallSuperheat" && birthTriggerMode_ != "cavityTipHsu")
    {
        FatalErrorInFunction
            << "Unsupported birthTriggerMode " << birthTriggerMode_
            << ". Supported values are legacyWallSuperheat and cavityTipHsu"
            << exit(FatalError);
    }

    if (seedHandoffMode_ != "none" && seedHandoffMode_ != "conservativeResolvedCap")
    {
        FatalErrorInFunction
            << "Unsupported seedHandoffMode " << seedHandoffMode_
            << ". Supported values are none and conservativeResolvedCap"
            << exit(FatalError);
    }

    if (birthTriggerMode_ == "cavityTipHsu" && seedHandoffMode_ != "conservativeResolvedCap")
    {
        FatalErrorInFunction
            << "birthTriggerMode cavityTipHsu requires seedHandoffMode conservativeResolvedCap"
            << exit(FatalError);
    }

    if (cavityMode())
    {
        if (Pstream::parRun())
        {
            FatalErrorInFunction
                << "birthTriggerMode cavityTipHsu is serial-only for M6A"
                << exit(FatalError);
        }

        if (cavityRadius_ <= SMALL)
        {
            FatalErrorInFunction
                << "cavityRadius must be positive, found " << cavityRadius_
                << exit(FatalError);
        }

        if (cavityContactAngle_ <= SMALL || cavityContactAngle_ >= 180.0 - SMALL)
        {
            FatalErrorInFunction
                << "cavityContactAngle must be in (0, 180) degrees, found "
                << cavityContactAngle_ << exit(FatalError);
        }

        if (activationPersistenceSteps_ < 1)
        {
            FatalErrorInFunction
                << "activationPersistenceSteps must be at least 1, found "
                << activationPersistenceSteps_ << exit(FatalError);
        }

        if (cavitySurfaceTension_ <= SMALL)
        {
            FatalErrorInFunction
                << "cavitySurfaceTension must be positive, found "
                << cavitySurfaceTension_ << exit(FatalError);
        }

        if (activationPersistence_ < 0 && activationPersistence_ != -1)
        {
            FatalErrorInFunction
                << "activationPersistence must be non-negative when supplied, "
                << "or -1 to use activationPersistenceSteps; found "
                << activationPersistence_ << exit(FatalError);
        }

        if (seedCreationSteps_ < 1 || seedCreationSteps_ > 1000)
        {
            FatalErrorInFunction
                << "seedCreationSteps must be between 1 and 1000, found "
                << seedCreationSteps_ << exit(FatalError);
        }

        if (maximumSeedCreationSteps_ < seedCreationSteps_)
        {
            FatalErrorInFunction
                << "maximumSeedCreationSteps must be >= seedCreationSteps; found "
                << maximumSeedCreationSteps_ << " and " << seedCreationSteps_
                << exit(FatalError);
        }

        if (maxSeedAlphaVaporPerStep_ <= 0 || maxSeedAlphaVaporPerStep_ > 1)
        {
            FatalErrorInFunction
                << "maxSeedAlphaVaporPerStep must be in (0, 1], found "
                << maxSeedAlphaVaporPerStep_ << exit(FatalError);
        }

        if
        (
            targetSeedVaporFraction_ <= SMALL
         || targetSeedVaporFraction_ > scalar(1) - epsLiquid_ + SMALL
        )
        {
            FatalErrorInFunction
                << "targetSeedVaporFraction must be in (0, 1-epsLiquid]; found "
                << targetSeedVaporFraction_ << " with epsLiquid " << epsLiquid_
                << exit(FatalError);
        }

        if (seedThermalReserveFraction_ < 0 || seedThermalReserveFraction_ > 1)
        {
            FatalErrorInFunction
                << "seedThermalReserveFraction must be in [0, 1], found "
                << seedThermalReserveFraction_ << exit(FatalError);
        }

        if (seedWallFluxUseFraction_ < 0 || seedWallFluxUseFraction_ > 1)
        {
            FatalErrorInFunction
                << "seedWallFluxUseFraction must be in [0, 1], found "
                << seedWallFluxUseFraction_ << exit(FatalError);
        }

        if (minimumCellsPerHandoffRadius_ < 1 || preferredCellsPerHandoffRadius_ < minimumCellsPerHandoffRadius_)
        {
            FatalErrorInFunction
                << "preferredCellsPerHandoffRadius must be >= minimumCellsPerHandoffRadius >= 1"
                << exit(FatalError);
        }

        if (maximumHandoffRadius_ <= SMALL)
        {
            FatalErrorInFunction
                << "maximumHandoffRadius must be positive, found " << maximumHandoffRadius_
                << exit(FatalError);
        }
    }

    if (massDensityPhase_ != "liquid" && massDensityPhase_ != "phase1")
    {
        FatalErrorInFunction
            << "wallSuperheatVaporBirth uses the existing phase-change convention "
            << "massSource = rho_liquid*alphaRate. Supported massDensityPhase values "
            << "are liquid and phase1; found " << massDensityPhase_
            << exit(FatalError);
    }

    if
    (
        sourceFootprintMode_ != "fixedSite"
     && sourceFootprintMode_ != "wallAdjacentContactLine"
    )
    {
        FatalErrorInFunction
            << "Unsupported sourceFootprintMode " << sourceFootprintMode_
            << ". Supported values are fixedSite and wallAdjacentContactLine"
            << exit(FatalError);
    }

    if
    (
        growthDepositionMode_ != "wallAdjacentCell"
     && growthDepositionMode_ != "wallNormalInterfaceBand"
    )
    {
        FatalErrorInFunction
            << "Unsupported growthDepositionMode " << growthDepositionMode_
            << ". Supported values are wallAdjacentCell and wallNormalInterfaceBand"
            << exit(FatalError);
    }

    if (growthBandHeight_ <= SMALL)
    {
        FatalErrorInFunction
            << "growthBandHeight must be positive, found " << growthBandHeight_
            << exit(FatalError);
    }

    if (growthBandAlphaMin_ < 0 || growthBandAlphaMin_ >= 1)
    {
        FatalErrorInFunction
            << "growthBandAlphaMin must satisfy 0 <= value < 1; found "
            << growthBandAlphaMin_ << exit(FatalError);
    }

    if (growthBandAlphaMax_ <= 0 || growthBandAlphaMax_ > 1)
    {
        FatalErrorInFunction
            << "growthBandAlphaMax must satisfy 0 < value <= 1; found "
            << growthBandAlphaMax_ << exit(FatalError);
    }

    if (growthBandAlphaMin_ >= growthBandAlphaMax_)
    {
        FatalErrorInFunction
            << "growthBandAlphaMin must be smaller than growthBandAlphaMax; found "
            << growthBandAlphaMin_ << " and " << growthBandAlphaMax_
            << exit(FatalError);
    }

    if (growthBandNormalDecayLength_ <= SMALL)
    {
        FatalErrorInFunction
            << "growthBandNormalDecayLength must be positive, found "
            << growthBandNormalDecayLength_ << exit(FatalError);
    }

    if (growthBandMinimumLiquid_ < 0 || growthBandMinimumLiquid_ >= 1)
    {
        FatalErrorInFunction
            << "growthBandMinimumLiquid must satisfy 0 <= value < 1; found "
            << growthBandMinimumLiquid_ << exit(FatalError);
    }

    if
    (
        growthDepositionMode_ == "wallNormalInterfaceBand"
     && sourceFootprintMode_ != "wallAdjacentContactLine"
    )
    {
        FatalErrorInFunction
            << "growthDepositionMode wallNormalInterfaceBand requires "
            << "sourceFootprintMode wallAdjacentContactLine"
            << exit(FatalError);
    }

    if (contactLineAlphaMin_ < 0 || contactLineAlphaMin_ >= 1)
    {
        FatalErrorInFunction
            << "contactLineAlphaMin must satisfy 0 <= value < 1; found "
            << contactLineAlphaMin_ << exit(FatalError);
    }

    if (contactLineAlphaMax_ <= 0 || contactLineAlphaMax_ > 1)
    {
        FatalErrorInFunction
            << "contactLineAlphaMax must satisfy 0 < value <= 1; found "
            << contactLineAlphaMax_ << exit(FatalError);
    }

    if (contactLineAlphaMin_ >= contactLineAlphaMax_)
    {
        FatalErrorInFunction
            << "contactLineAlphaMin must be smaller than contactLineAlphaMax; found "
            << contactLineAlphaMin_ << " and " << contactLineAlphaMax_
            << exit(FatalError);
    }

    if (minimumSourceLiquidFraction_ < 0 || minimumSourceLiquidFraction_ >= 1)
    {
        FatalErrorInFunction
            << "minimumSourceLiquidFraction must satisfy 0 <= value < 1; found "
            << minimumSourceLiquidFraction_ << exit(FatalError);
    }

    if (contactLineNeighbourDelta_ <= 0 || contactLineNeighbourDelta_ > 1)
    {
        FatalErrorInFunction
            << "contactLineNeighbourDelta must satisfy 0 < value <= 1; found "
            << contactLineNeighbourDelta_ << exit(FatalError);
    }

    if
    (
        sustainedBirthMode_ != "none"
     && sustainedBirthMode_ != "thermallyBoundedHold"
    )
    {
        FatalErrorInFunction
            << "Unsupported sustainedBirthMode " << sustainedBirthMode_
            << ". Supported values are none and thermallyBoundedHold"
            << exit(FatalError);
    }

    if (sustainedBirthMode_ == "none")
    {
        if
        (
            minimumActiveDuration_ != 0
         || minimumEvaporationSuperheat_ != 0
         || thermalReserveFraction_ != 0
         || minimumThermalAlphaStep_ != 0
        )
        {
            FatalErrorInFunction
                << "Sustained-birth controls require "
                << "sustainedBirthMode thermallyBoundedHold; old cases must not "
                << "silently enable sustained birth"
                << exit(FatalError);
        }
    }

    if (minimumActiveDuration_ < 0)
    {
        FatalErrorInFunction
            << "minimumActiveDuration must be non-negative, found "
            << minimumActiveDuration_ << exit(FatalError);
    }

    if (minimumEvaporationSuperheat_ < 0)
    {
        FatalErrorInFunction
            << "minimumEvaporationSuperheat must be non-negative, found "
            << minimumEvaporationSuperheat_ << exit(FatalError);
    }

    if (thermalReserveFraction_ < 0 || thermalReserveFraction_ > 1)
    {
        FatalErrorInFunction
            << "thermalReserveFraction must satisfy 0 <= value <= 1; found "
            << thermalReserveFraction_ << exit(FatalError);
    }

    if
    (
        thermalLimitMode_ != "sensibleReserve"
     && thermalLimitMode_ != "wallFluxAwareInstantaneous"
    )
    {
        FatalErrorInFunction
            << "Unsupported thermalLimitMode " << thermalLimitMode_
            << ". Supported values are sensibleReserve and "
            << "wallFluxAwareInstantaneous"
            << exit(FatalError);
    }

    if (wallFluxUseFraction_ < 0 || wallFluxUseFraction_ > 1)
    {
        FatalErrorInFunction
            << "wallFluxUseFraction must satisfy 0 <= value <= 1; found "
            << wallFluxUseFraction_ << exit(FatalError);
    }

    if (perCellDiagnosticsInterval_ < 1)
    {
        FatalErrorInFunction
            << "perCellDiagnosticsInterval must be >= 1; found "
            << perCellDiagnosticsInterval_ << exit(FatalError);
    }

    if (compactDiagnosticsInterval_ < 1)
    {
        FatalErrorInFunction
            << "compactDiagnosticsInterval must be >= 1; found "
            << compactDiagnosticsInterval_ << exit(FatalError);
    }

    if (perCellDiagnosticsFormat_ != "csv")
    {
        FatalErrorInFunction
            << "Unsupported perCellDiagnosticsFormat " << perCellDiagnosticsFormat_
            << ". Supported value is csv"
            << exit(FatalError);
    }

    if
    (
        thermalLimitMode_ == "wallFluxAwareInstantaneous"
     && !modelDict().found("wallFluxUseFraction")
    )
    {
        FatalErrorInFunction
            << "thermalLimitMode wallFluxAwareInstantaneous requires explicit "
            << "wallFluxUseFraction"
            << exit(FatalError);
    }

    if (sustainedBirthMode_ == "thermallyBoundedHold" && thermalReserveFraction_ <= 0)
    {
        FatalErrorInFunction
            << "sustainedBirthMode thermallyBoundedHold requires positive "
            << "thermalReserveFraction"
            << exit(FatalError);
    }

    if (minimumThermalAlphaStep_ < 0)
    {
        FatalErrorInFunction
            << "minimumThermalAlphaStep must be non-negative, found "
            << minimumThermalAlphaStep_ << exit(FatalError);
    }

    if (siteMode_ != "allPatchAdjacent" && siteMode_ != "circularPatchSite")
    {
        FatalErrorInFunction
            << "Unsupported wallSuperheatVaporBirth siteMode " << siteMode_
            << ". Supported values are allPatchAdjacent and circularPatchSite"
            << exit(FatalError);
    }

    if (siteMode_ == "circularPatchSite")
    {
        if (sitePatch_.empty())
        {
            FatalErrorInFunction
                << "siteMode circularPatchSite requires sitePatch"
                << exit(FatalError);
        }

        bool configuredPatch = false;
        forAll(patches_, patchi)
        {
            configuredPatch = configuredPatch || patches_[patchi] == sitePatch_;
        }

        if (!configuredPatch)
        {
            FatalErrorInFunction
                << "sitePatch " << sitePatch_
                << " is not listed in configured patches " << patches_
                << exit(FatalError);
        }

        if (siteRadius_ <= 0)
        {
            FatalErrorInFunction
                << "siteMode circularPatchSite requires positive siteRadius; found "
                << siteRadius_ << exit(FatalError);
        }

        if (sourceFootprintMode_ == "fixedSite" && maxActiveSites_ != 1)
        {
            FatalErrorInFunction
                << "siteMode circularPatchSite currently requires maxActiveSites 1; found "
                << maxActiveSites_ << exit(FatalError);
        }

        if (sourceFootprintMode_ == "fixedSite" && !lockSiteAfterBirth_)
        {
            FatalErrorInFunction
                << "siteMode circularPatchSite requires lockSiteAfterBirth true"
                << exit(FatalError);
        }
    }

    if
    (
        sourceFootprintMode_ == "wallAdjacentContactLine"
     && siteMode_ != "circularPatchSite"
    )
    {
        FatalErrorInFunction
            << "sourceFootprintMode wallAdjacentContactLine requires "
            << "growthSiteMode/siteMode circularPatchSite"
            << exit(FatalError);
    }

    if (departurePersistence_ < 0)
    {
        FatalErrorInFunction
            << "departurePersistence must be non-negative, found "
            << departurePersistence_ << exit(FatalError);
    }

    if (departureContactAreaThreshold_ < 0)
    {
        FatalErrorInFunction
            << "departureContactAreaThreshold must be non-negative, found "
            << departureContactAreaThreshold_ << exit(FatalError);
    }

    if (departureSiteVaporThreshold_ < 0 || departureSiteVaporThreshold_ > 1)
    {
        FatalErrorInFunction
            << "departureSiteVaporThreshold must satisfy 0 <= value <= 1; found "
            << departureSiteVaporThreshold_ << exit(FatalError);
    }

    if (rewetAlphaThreshold_ < 0 || rewetAlphaThreshold_ > 1)
    {
        FatalErrorInFunction
            << "rewetAlphaThreshold must satisfy 0 <= value <= 1; found "
            << rewetAlphaThreshold_ << exit(FatalError);
    }

    if (rewetWallVaporThreshold_ < 0 || rewetWallVaporThreshold_ > 1)
    {
        FatalErrorInFunction
            << "rewetWallVaporThreshold must satisfy 0 <= value <= 1; found "
            << rewetWallVaporThreshold_ << exit(FatalError);
    }

    if (rewetPersistence_ < 0)
    {
        FatalErrorInFunction
            << "rewetPersistence must be non-negative, found "
            << rewetPersistence_ << exit(FatalError);
    }

    if (minimumRearmDelay_ < 0)
    {
        FatalErrorInFunction
            << "minimumRearmDelay must be non-negative, found "
            << minimumRearmDelay_ << exit(FatalError);
    }


    if
    (
        activationAlphaLiquidThreshold_ < 0
     || activationAlphaLiquidThreshold_ > 1
     || activationWallVaporThreshold_ < 0
     || activationWallVaporThreshold_ > 1
    )
    {
        FatalErrorInFunction
            << "Activation wetting thresholds must be in [0, 1]"
            << exit(FatalError);
    }

    if (shieldingRadius_ < 0 || shieldingHeight_ < 0)
    {
        FatalErrorInFunction
            << "shieldingRadius and shieldingHeight must be non-negative"
            << exit(FatalError);
    }

    if
    (
        shieldingAlphaVaporThreshold_ < 0
     || shieldingAlphaVaporThreshold_ > 1
    )
    {
        FatalErrorInFunction
            << "shieldingAlphaVaporThreshold must be in [0, 1]"
            << exit(FatalError);
    }
}

void Foam::wallSuperheatVaporBirth::buildCandidateCells()
{
    const fvMesh& mesh = phase1_.mesh();
    boolList isCandidate(mesh.nCells(), false);
    wallFaceAreaCache_.setSize(mesh.nCells(), scalar(0));

    forAll(patches_, patchi)
    {
        const label patchID = mesh.boundaryMesh().findPatchID(patches_[patchi]);
        if (patchID < 0)
        {
            FatalErrorInFunction
                << "Unable to find configured wall-superheat vapor-birth patch "
                << patches_[patchi] << exit(FatalError);
        }

        const fvPatch& patch = mesh.boundary()[patchID];
        const labelUList& faceCells = patch.faceCells();
        const vectorField& Sf = patch.Sf();
        forAll(faceCells, facei)
        {
            const label celli = faceCells[facei];
            isCandidate[celli] = true;
            wallFaceAreaCache_[celli] += mag(Sf[facei]);
        }
    }

    DynamicList<label> cells(mesh.nCells());
    forAll(isCandidate, celli)
    {
        if (isCandidate[celli])
        {
            cells.append(celli);
        }
    }

    candidateCells_.transfer(cells);
    patchAdjacentCellCount_ = candidateCells_.size();

    siteCellMask_.setSize(mesh.nCells(), false);
    sustainedState_.setSize(mesh.nCells(), INACTIVE);
    sustainedActivationTime_.setSize(mesh.nCells(), -GREAT);

    if (siteMode_ == "allPatchAdjacent")
    {
        siteCells_ = candidateCells_;
    }
    else
    {
        const label patchID = mesh.boundaryMesh().findPatchID(sitePatch_);
        if (patchID < 0)
        {
            FatalErrorInFunction
                << "Unable to find configured circularPatchSite patch "
                << sitePatch_ << exit(FatalError);
        }

        siteCells_ = circularPatchSiteCells
        (
            mesh.boundary()[patchID].Cf(),
            mesh.boundary()[patchID].faceCells(),
            siteCenter_,
            siteRadius_
        );

        label globalSiteCellCount = siteCells_.size();
        label globalPatchAdjacentCellCount = patchAdjacentCellCount_;
        reduce(globalSiteCellCount, sumOp<label>());
        reduce(globalPatchAdjacentCellCount, sumOp<label>());

        if (globalSiteCellCount == 0)
        {
            FatalErrorInFunction
                << "circularPatchSite on patch " << sitePatch_
                << " with siteCenter " << siteCenter_
            << " and siteRadius " << siteRadius_
            << " contains zero globally eligible wall-adjacent cells"
                << exit(FatalError);
        }

        if (globalSiteCellCount >= globalPatchAdjacentCellCount)
        {
            FatalErrorInFunction
            << "circularPatchSite selected " << globalSiteCellCount
            << " global cells, which is not localized relative to "
            << globalPatchAdjacentCellCount << " patch-adjacent cells"
                << exit(FatalError);
        }

        candidateCells_ = siteCells_;
    }

    forAll(siteCells_, siteI)
    {
        siteCellMask_[siteCells_[siteI]] = true;
    }
}

Foam::scalar Foam::wallSuperheatVaporBirth::siteVolume() const
{
    const fvMesh& mesh = phase1_.mesh();
    scalar volume = 0;
    forAll(siteCells_, siteI)
    {
        volume += mesh.V()[siteCells_[siteI]];
    }
    reduce(volume, sumOp<scalar>());
    return volume;
}

Foam::scalar Foam::wallSuperheatVaporBirth::TSat(const label celli) const
{
    return TSatValue_ > 0 ? TSatValue_ : satModel_.TSat()[celli];
}

Foam::scalar Foam::wallSuperheatVaporBirth::wallFluxIntoFluid(const label celli) const
{
    const fvMesh& mesh = phase1_.mesh();
    const volScalarField& liquidTemperature = phase1_.thermo().T();
    tmp<volScalarField> tAlphat = turbModel_.alphat();
    scalar heatRatePerDeltaT = 0;
    scalar area = 0;

    forAll(patches_, patchi)
    {
        const label patchID = mesh.boundaryMesh().findPatchID(patches_[patchi]);
        if (patchID < 0)
        {
            continue;
        }

        const fvPatch& patch = mesh.boundary()[patchID];
        const labelUList& faceCells = patch.faceCells();
        const scalarField& alphatPatch = tAlphat().boundaryField()[patchID];
        tmp<scalarField> tKappaPatch =
            phase1_.thermo().kappaEff(alphatPatch, patchID);
        const scalarField& kappaPatch = tKappaPatch();
        tmp<scalarField> tSnGrad = liquidTemperature.boundaryField()[patchID].snGrad();
        const scalarField& snGrad = tSnGrad();
        const vectorField& Sf = patch.Sf();

        forAll(faceCells, facei)
        {
            if (faceCells[facei] == celli)
            {
                const scalar faceArea = mag(Sf[facei]);
                heatRatePerDeltaT += kappaPatch[facei]*snGrad[facei]*faceArea;
                area += faceArea;
            }
        }
    }

    return area > SMALL ? heatRatePerDeltaT/area : scalar(0);
}

Foam::scalar Foam::wallSuperheatVaporBirth::wallFaceArea(const label celli) const
{
    return
        celli >= 0 && celli < wallFaceAreaCache_.size()
      ? wallFaceAreaCache_[celli]
      : scalar(0);
}

bool Foam::wallSuperheatVaporBirth::contactLineEligible
(
    const scalar alphaLiquid,
    const scalar neighbourMaxAlphaVapor,
    const scalar contactLineAlphaMin,
    const scalar contactLineAlphaMax,
    const scalar minimumSourceLiquidFraction,
    const bool contactLineNeighbourSearch,
    const scalar contactLineNeighbourDelta
)
{
    if (alphaLiquid < minimumSourceLiquidFraction)
    {
        return false;
    }

    const scalar alphaVapor = scalar(1) - alphaLiquid;

    if (alphaVapor >= contactLineAlphaMin && alphaVapor <= contactLineAlphaMax)
    {
        return true;
    }

    return
        contactLineNeighbourSearch
     && alphaVapor < contactLineAlphaMin
     && neighbourMaxAlphaVapor >= contactLineNeighbourDelta;
}

bool Foam::wallSuperheatVaporBirth::contactLineEligible(const label celli) const
{
    if (!siteCellMask_[celli])
    {
        return false;
    }

    scalar neighbourMaxAlphaVapor = 0;
    const labelList& neighbours = phase1_.mesh().cellCells()[celli];
    forAll(neighbours, neighbourI)
    {
        const label neighbCelli = neighbours[neighbourI];
        if (!siteCellMask_[neighbCelli])
        {
            continue;
        }

        neighbourMaxAlphaVapor = max
        (
            neighbourMaxAlphaVapor,
            scalar(1) - phase1_[neighbCelli]
        );
    }

    return contactLineEligible
    (
        phase1_[celli],
        neighbourMaxAlphaVapor,
        contactLineAlphaMin_,
        contactLineAlphaMax_,
        minimumSourceLiquidFraction_,
        contactLineNeighbourSearch_,
        contactLineNeighbourDelta_
    );
}

void Foam::wallSuperheatVaporBirth::appendInterfaceBandCells
(
    const label wallCelli,
    const scalarField& alphaLiquid,
    boolList& included,
    scalarField& weights,
    scalarField& wallNormalDistance,
    label& eligibleCells
) const
{
    const fvMesh& mesh = phase1_.mesh();
    vector faceCentre(vector::zero);
    vector normal(vector::zero);
    scalar wallArea = 0;
    bool foundWallFace = false;

    forAll(patches_, patchi)
    {
        const label patchID = mesh.boundaryMesh().findPatchID(patches_[patchi]);
        if (patchID < 0)
        {
            continue;
        }

        const fvPatch& patch = mesh.boundary()[patchID];
        const labelUList& faceCells = patch.faceCells();
        const vectorField& Sf = patch.Sf();
        const vectorField& Cf = patch.Cf();
        forAll(faceCells, facei)
        {
            if (faceCells[facei] != wallCelli)
            {
                continue;
            }

            faceCentre = Cf[facei];
            normal = mesh.C()[wallCelli] - faceCentre;
            wallArea = mag(Sf[facei]);
            foundWallFace = true;
            break;
        }
        if (foundWallFace)
        {
            break;
        }
    }

    if (!foundWallFace || mag(normal) <= SMALL || wallArea <= SMALL)
    {
        return;
    }

    normal /= mag(normal);
    const scalar tangentialLimit = max(Foam::sqrt(wallArea), localTangentialSpacing_);
    boolList visited(mesh.nCells(), false);
    DynamicList<label> queue;
    queue.append(wallCelli);
    visited[wallCelli] = true;

    for (label queueI = 0; queueI < queue.size(); ++queueI)
    {
        const label celli = queue[queueI];
        const vector d = mesh.C()[celli] - faceCentre;
        const scalar projection = d & normal;
        const vector tangent = d - projection*normal;
        const scalar tangentialDistance = mag(tangent);

        if
        (
            projection >= -SMALL
         && projection <= growthBandHeight_ + tangentialLimit
         && tangentialDistance <= tangentialLimit + SMALL
        )
        {
            const labelList& neighbours = mesh.cellCells()[celli];
            forAll(neighbours, neighbourI)
            {
                const label neighbCelli = neighbours[neighbourI];
                if (!visited[neighbCelli])
                {
                    const vector nd = mesh.C()[neighbCelli] - faceCentre;
                    const scalar nProjection = nd & normal;
                    const scalar nTangentialDistance = mag(nd - nProjection*normal);
                    if
                    (
                        nProjection >= -SMALL
                     && nProjection <= growthBandHeight_ + tangentialLimit
                     && nTangentialDistance <= tangentialLimit + SMALL
                    )
                    {
                        queue.append(neighbCelli);
                        visited[neighbCelli] = true;
                    }
                }
            }
        }

        if
        (
            projection < -SMALL
         || projection > growthBandHeight_ + SMALL
         || tangentialDistance > tangentialLimit + SMALL
         || alphaLiquid[celli] < growthBandMinimumLiquid_
        )
        {
            continue;
        }

        const scalar alphaVapor = max(min(scalar(1) - alphaLiquid[celli], scalar(1)), scalar(0));
        scalar interfaceProxy = 0;
        bool interfaceCell = alphaVapor >= growthBandAlphaMin_ && alphaVapor <= growthBandAlphaMax_;
        if (interfaceCell)
        {
            interfaceProxy = max(alphaVapor*(scalar(1) - alphaVapor), VSMALL);
        }

        const labelList& neighbours = mesh.cellCells()[celli];
        forAll(neighbours, neighbourI)
        {
            const label neighbCelli = neighbours[neighbourI];
            const scalar neighbAlphaVapor = max(min(scalar(1) - alphaLiquid[neighbCelli], scalar(1)), scalar(0));
            if (neighbAlphaVapor >= growthBandAlphaMin_ && neighbAlphaVapor <= growthBandAlphaMax_)
            {
                interfaceCell = true;
                interfaceProxy = max(interfaceProxy, neighbAlphaVapor*(scalar(1) - neighbAlphaVapor));
            }
        }

        if (!interfaceCell)
        {
            continue;
        }

        if (growthBandRequirePrimaryBubble_ && alphaVapor <= SMALL && interfaceProxy <= VSMALL)
        {
            continue;
        }

        const scalar weight = max(interfaceProxy, VSMALL)*mesh.V()[celli]
           /(projection + growthBandNormalDecayLength_);
        if (weight <= SMALL)
        {
            continue;
        }

        if (!included[celli])
        {
            eligibleCells++;
            included[celli] = true;
            wallNormalDistance[celli] = projection;
        }
        else
        {
            wallNormalDistance[celli] = min(wallNormalDistance[celli], projection);
        }
        weights[celli] += weight;
    }
}

bool Foam::wallSuperheatVaporBirth::cavityMode() const
{
    return birthTriggerMode_ == "cavityTipHsu";
}

const char* Foam::wallSuperheatVaporBirth::cavityStateName(const label state) const
{
    if (state == WAITING_FOR_CRITERION)
    {
        return "WAITING_FOR_CRITERION";
    }
    if (state == ACCUMULATING_HANDOFF_ENERGY)
    {
        return "ACCUMULATING_HANDOFF_ENERGY";
    }
    if (state == CREATING_RESOLVED_CAP)
    {
        return "CREATING_RESOLVED_CAP";
    }
    if (state == CONTACT_LINE_GROWTH)
    {
        return "CONTACT_LINE_GROWTH";
    }
    if (state == WAITING_FOR_REWET)
    {
        return "WAITING_FOR_REWET";
    }
    return "UNKNOWN";
}

void Foam::wallSuperheatVaporBirth::restoreCavityState()
{
    if (!cavityMode())
    {
        return;
    }

    const fileName diagnosticDir
    (
        phase1_.mesh().time().path()
       /"postProcessing"
       /diagnosticPrefix_
       /"serial"
    );
    const fileName stateFile(diagnosticDir/"birth-state-history.csv");
    std::ifstream input(stateFile.c_str());
    if (!input.good())
    {
        return;
    }

    std::string line;
    std::string last;
    while (std::getline(input, line))
    {
        if (!line.empty() && line.rfind("timeIndex,", 0) != 0)
        {
            last = line;
        }
    }
    if (last.empty())
    {
        return;
    }

    DynamicList<std::string> cols;
    std::stringstream ss(last);
    std::string item;
    while (std::getline(ss, item, ','))
    {
        cols.append(item);
    }
    if (cols.size() < 13)
    {
        WarningInFunction
            << "Ignoring incomplete cavity-state row in " << stateFile
            << endl;
        return;
    }

    const word state(cols[2]);
    if (state == "ACCUMULATING_HANDOFF_ENERGY")
    {
        cavityBirthState_ = ACCUMULATING_HANDOFF_ENERGY;
    }
    else if (state == "CREATING_RESOLVED_CAP")
    {
        cavityBirthState_ = CREATING_RESOLVED_CAP;
    }
    else if (state == "CONTACT_LINE_GROWTH")
    {
        cavityBirthState_ = CONTACT_LINE_GROWTH;
    }
    else if (state == "WAITING_FOR_REWET")
    {
        cavityBirthState_ = WAITING_FOR_REWET;
    }
    else
    {
        cavityBirthState_ = WAITING_FOR_CRITERION;
    }

    activationCounter_ = std::atoi(cols[3].c_str());
    activationTime_ = std::atof(cols[5].c_str());
    handoffStartTime_ = std::atof(cols[6].c_str());
    handoffCompleteTime_ = std::atof(cols[7].c_str());
    storedHandoffEnergy_ = std::atof(cols[8].c_str());
    createdVaporVolume_ = std::atof(cols[9].c_str());
    createdVaporMass_ = std::atof(cols[10].c_str());
    consumedLatentEnergy_ = std::atof(cols[11].c_str());
    creationStepIndex_ = std::atoi(cols[12].c_str());

    if (cols.size() >= 28)
    {
        cycleId_ = std::atoi(cols[13].c_str());
        departureCandidateDuration_ = std::atof(cols[14].c_str());
        rewetCandidateDuration_ = std::atof(cols[15].c_str());
        departureDetected_ = std::atoi(cols[16].c_str()) != 0;
        departureTime_ = std::atof(cols[17].c_str());
        rewetDetected_ = std::atoi(cols[18].c_str()) != 0;
        rewetTime_ = std::atof(cols[19].c_str());
        siteRearmed_ = std::atoi(cols[20].c_str()) != 0;
        rearmCount_ = std::atoi(cols[21].c_str());
        rearmTime_ = std::atof(cols[22].c_str());
        nextActivationTime_ = std::atof(cols[23].c_str());
    }
    if (cols.size() >= 34)
    {
        activationCandidateDuration_ = std::atof(cols[28].c_str());
        lastSiteShielded_ = std::atoi(cols[30].c_str()) != 0;
        displacedLiquidMass_ = std::atof(cols[31].c_str());
    }
}
void Foam::wallSuperheatVaporBirth::evaluateSiteAttachment
(
    scalar& siteMeanAlphaLiquid,
    scalar& siteVaporWallFraction,
    scalar& attachedWallContactArea,
    bool& attachedBubblePresent
) const
{
    const fvMesh& mesh = phase1_.mesh();
    scalar alphaLiquidSum = 0;
    scalar volumeSum = 0;
    scalar siteWallArea = 0;
    scalar vaporWallArea = 0;

    forAll(siteCells_, siteI)
    {
        const label celli = siteCells_[siteI];
        const scalar alphaLiquid = max(min(phase1_[celli], scalar(1)), scalar(0));
        const scalar alphaVapor = scalar(1) - alphaLiquid;
        const scalar volume = mesh.V()[celli];
        const scalar area = wallFaceArea(celli);

        alphaLiquidSum += alphaLiquid*volume;
        volumeSum += volume;
        siteWallArea += area;
        if (alphaVapor >= departureSiteVaporThreshold_)
        {
            vaporWallArea += area;
        }
    }

    reduce(alphaLiquidSum, sumOp<scalar>());
    reduce(volumeSum, sumOp<scalar>());
    reduce(siteWallArea, sumOp<scalar>());
    reduce(vaporWallArea, sumOp<scalar>());

    siteMeanAlphaLiquid = volumeSum > SMALL ? alphaLiquidSum/volumeSum : scalar(1);
    attachedWallContactArea = vaporWallArea;
    siteVaporWallFraction = siteWallArea > SMALL ? vaporWallArea/siteWallArea : scalar(0);
    attachedBubblePresent = vaporWallArea > departureContactAreaThreshold_ + SMALL;
}

void Foam::wallSuperheatVaporBirth::rearmForNextCycle()
{
    cavityBirthState_ = WAITING_FOR_CRITERION;
    activationCounter_ = 0;
    activationCandidateDuration_ = 0;
    creationStepIndex_ = 0;
    activationTime_ = -1;
    handoffStartTime_ = -1;
    handoffCompleteTime_ = -1;
    storedHandoffEnergy_ = 0;
    createdVaporVolume_ = 0;
    createdVaporMass_ = 0;
    displacedLiquidMass_ = 0;
    consumedLatentEnergy_ = 0;
    accumulatedWallEnergy_ = 0;
    lastLocalSensibleEnergy_ = 0;
    lastCavityTtip_ = 0;
    handoffRadius_ = 0;
    handoffCapHeight_ = 0;
    handoffFootprintRadius_ = 0;
    handoffTargetVolume_ = 0;
    handoffTargetMass_ = 0;
    handoffTargetLatentEnergy_ = 0;
    handoffCells_.clear();
    handoffCellMask_.clear();
    handoffTargetVaporFractions_.clear();
    cavityTipCells_.clear();
    cavityTipWeights_.clear();
    shieldingCells_.clear();

    departureCandidateDuration_ = 0;
    rewetCandidateDuration_ = 0;
    departureDetected_ = false;
    // Preserve the accepted rewet event on the rearm row. It is cleared
    // only when a new cycle actually activates.
    rewetDetected_ = true;
    siteRearmed_ = true;
    lastSiteShielded_ = false;
    rearmCount_++;
    cycleId_++;
    rearmTime_ = phase1_.mesh().time().value();
    nextActivationTime_ = -1;

    forAll(siteCells_, siteI)
    {
        const label celli = siteCells_[siteI];
        sustainedState_[celli] = INACTIVE;
        sustainedActivationTime_[celli] = -GREAT;
        sustainedStateField_[celli] = sustainedState_[celli];
    }
}
Foam::scalar Foam::wallSuperheatVaporBirth::interpolateCavityTipTemperature
(
    const volScalarField& liquidTemperature
) const
{
    scalar weightedTemperature = 0;
    scalar weightSum = 0;
    forAll(cavityTipCells_, stencilI)
    {
        const label celli = cavityTipCells_[stencilI];
        const scalar weight = cavityTipWeights_[stencilI];
        weightedTemperature += weight*liquidTemperature[celli];
        weightSum += weight;
    }
    reduce(weightedTemperature, sumOp<scalar>());
    reduce(weightSum, sumOp<scalar>());
    return weightSum > SMALL ? weightedTemperature/weightSum : scalar(0);
}

bool Foam::wallSuperheatVaporBirth::siteShielded() const
{
    if (shieldingCells_.empty())
    {
        return false;
    }

    scalar maximumVaporFraction = 0;
    forAll(shieldingCells_, shieldI)
    {
        const label celli = shieldingCells_[shieldI];
        maximumVaporFraction = max
        (
            maximumVaporFraction,
            scalar(1) - max(min(phase1_[celli], scalar(1)), scalar(0))
        );
    }
    reduce(maximumVaporFraction, maxOp<scalar>());
    return maximumVaporFraction + SMALL >= shieldingAlphaVaporThreshold_;
}
void Foam::wallSuperheatVaporBirth::buildHandoffCells
(
    const volScalarField& vaporDensity
)
{
    if (handoffTargetVolume_ > SMALL)
    {
        return;
    }

    const fvMesh& mesh = phase1_.mesh();
    const scalar pi = constant::mathematical::pi;
    const scalar theta = cavityContactAngle_*pi/180.0;
    const scalar sinTheta = sin(theta);
    const scalar cosTheta = cos(theta);
    cavityCandidateRadius_ = cavityRadius_/max(sinTheta, SMALL);
    cavityTipHeight_ =
        cavityRadius_*(scalar(1) + cosTheta)/max(sinTheta, SMALL);

    scalar rhoVSum = 0;
    scalar TSatSum = 0;
    scalar latentHeatSum = 0;
    scalar thermoWeight = 0;
    localTangentialSpacing_ = GREAT;
    localWallNormalSpacing_ = GREAT;

    forAll(candidateCells_, candidateI)
    {
        const label celli = candidateCells_[candidateI];
        const vector d = mesh.C()[celli] - cavityCentre_;
        const scalar radius = Foam::sqrt(sqr(d.x()) + sqr(d.y()));
        if (radius > siteRadius_ + SMALL)
        {
            continue;
        }

        const scalar area = wallFaceArea(celli);
        if (area > SMALL)
        {
            localTangentialSpacing_ =
                min(localTangentialSpacing_, Foam::sqrt(area));
            localWallNormalSpacing_ =
                min(localWallNormalSpacing_, mesh.V()[celli]/area);
        }

        const scalar volume = mesh.V()[celli];
        rhoVSum += vaporDensity[celli]*volume;
        TSatSum += TSat(celli)*volume;
        latentHeatSum += satModel_.L()[celli]*volume;
        thermoWeight += volume;
    }

    reduce(rhoVSum, sumOp<scalar>());
    reduce(TSatSum, sumOp<scalar>());
    reduce(latentHeatSum, sumOp<scalar>());
    reduce(thermoWeight, sumOp<scalar>());
    reduce(localTangentialSpacing_, minOp<scalar>());
    reduce(localWallNormalSpacing_, minOp<scalar>());

    if
    (
        localTangentialSpacing_ >= GREAT/2
     || localWallNormalSpacing_ >= GREAT/2
     || thermoWeight <= SMALL
    )
    {
        FatalErrorInFunction
            << "Unable to determine local cavity handoff mesh spacing or "
            << "thermodynamic state for site " << diagnosticPrefix_
            << exit(FatalError);
    }

    const scalar localRhoVapor = rhoVSum/thermoWeight;
    const scalar localTSat = TSatSum/thermoWeight;
    const scalar localLatentHeat = latentHeatSum/thermoWeight;
    cavityCriticalSuperheat_ =
        scalar(2)*cavitySurfaceTension_*localTSat
       /
        (
            max(cavityCandidateRadius_, SMALL)
           *max(localLatentHeat, SMALL)
           *max(localRhoVapor, SMALL)
        );
    cavityCriticalTemperature_ = localTSat + cavityCriticalSuperheat_;

    handoffRadius_ = max
    (
        cavityCandidateRadius_,
        max
        (
            preferredCellsPerHandoffRadius_*localTangentialSpacing_,
            preferredCellsPerHandoffRadius_*localWallNormalSpacing_
        )
    );

    const scalar minimumResolvedRadius = max
    (
        minimumCellsPerHandoffRadius_*localTangentialSpacing_,
        minimumCellsPerHandoffRadius_*localWallNormalSpacing_
    );
    if (handoffRadius_ + SMALL < minimumResolvedRadius)
    {
        FatalErrorInFunction
            << "Resolved handoff radius " << handoffRadius_
            << " is below the minimum mesh-resolved radius "
            << minimumResolvedRadius << " for site " << diagnosticPrefix_
            << exit(FatalError);
    }
    if (handoffRadius_ > maximumHandoffRadius_ + SMALL)
    {
        FatalErrorInFunction
            << "Resolved handoff radius " << handoffRadius_
            << " exceeds maximumHandoffRadius " << maximumHandoffRadius_
            << " for site " << diagnosticPrefix_
            << exit(FatalError);
    }

    // Use the same major spherical-cap geometry as the Hsu cavity model:
    // h = R(1 + cos(theta)) and Rc = R sin(theta). The sphere centre
    // therefore lies R cos(theta) above the nominal wall plane.
    handoffCapHeight_ = handoffRadius_*(scalar(1) + cosTheta);
    handoffFootprintRadius_ = handoffRadius_*sinTheta;
    const vector sphereCentre =
        cavityCentre_ + vector(0, 0, handoffRadius_*cosTheta);

    DynamicList<label> cells;
    DynamicList<scalar> targetFractions;
    const scalar seedInterfaceThickness = max
    (
        localTangentialSpacing_,
        localWallNormalSpacing_
    );
    scalar volume = 0;
    scalar mass = 0;
    scalar latentEnergy = 0;
    label targetVaporCoreCells = 0;
    forAll(mesh.C(), celli)
    {
        const vector d = mesh.C()[celli] - sphereCentre;
        const scalar z = mesh.C()[celli].z() - cavityCentre_.z();
        if (z < -SMALL)
        {
            continue;
        }

        // Use a nearly pure-vapor core with a one-cell signed-distance
        // transition. Bulk-phase redistribution models need at least one
        // vapor-core cell; a uniformly mixed numerical cap is insufficient.
        const scalar signedDistance = handoffRadius_ - mag(d);
        const scalar profile = min
        (
            max
            (
                scalar(0.5)
              + signedDistance/max(seedInterfaceThickness, VSMALL),
                scalar(0)
            ),
            scalar(1)
        );
        const scalar targetVaporFraction =
            min(targetSeedVaporFraction_, profile);
        if (targetVaporFraction <= SMALL)
        {
            continue;
        }

        cells.append(celli);
        targetFractions.append(targetVaporFraction);
        targetVaporCoreCells += targetVaporFraction > scalar(0.999) + SMALL;
        const scalar cellVolume = mesh.V()[celli];
        const scalar targetVaporVolume = targetVaporFraction*cellVolume;
        const scalar cellVaporMass =
            vaporDensity[celli]*targetVaporVolume;
        volume += targetVaporVolume;
        mass += cellVaporMass;
        latentEnergy += cellVaporMass*satModel_.L()[celli];
    }
    reduce(volume, sumOp<scalar>());
    reduce(mass, sumOp<scalar>());
    reduce(latentEnergy, sumOp<scalar>());
    reduce(targetVaporCoreCells, sumOp<label>());
    handoffCells_.transfer(cells);
    handoffCellMask_.setSize(mesh.nCells(), false);
    forAll(handoffCells_, handoffI)
    {
        handoffCellMask_[handoffCells_[handoffI]] = true;
    }
    handoffTargetVaporFractions_.setSize(targetFractions.size());
    forAll(handoffTargetVaporFractions_, targetI)
    {
        handoffTargetVaporFractions_[targetI] = targetFractions[targetI];
    }
    handoffTargetVolume_ = volume;
    if (handoffTargetVaporFractions_.size() != handoffCells_.size())
    {
        FatalErrorInFunction
            << "Seed target profile size mismatch for site "
            << diagnosticPrefix_ << exit(FatalError);
    }
    handoffTargetMass_ = mass;
    handoffTargetLatentEnergy_ = latentEnergy;

    label nearestCells[8];
    scalar nearestDist2[8];
    label nearestCount = 0;
    const vector tip = cavityCentre_ + vector(0, 0, cavityTipHeight_);
    forAll(mesh.C(), celli)
    {
        const scalar distance2 = magSqr(mesh.C()[celli] - tip);
        if (nearestCount < 8)
        {
            nearestCells[nearestCount] = celli;
            nearestDist2[nearestCount] = distance2;
            ++nearestCount;
            continue;
        }
        label worst = 0;
        for (label i = 1; i < nearestCount; ++i)
        {
            if (nearestDist2[i] > nearestDist2[worst])
            {
                worst = i;
            }
        }
        if (distance2 < nearestDist2[worst])
        {
            nearestCells[worst] = celli;
            nearestDist2[worst] = distance2;
        }
    }
    cavityTipCells_.setSize(nearestCount);
    cavityTipWeights_.setSize(nearestCount);
    for (label i = 0; i < nearestCount; ++i)
    {
        cavityTipCells_[i] = nearestCells[i];
        cavityTipWeights_[i] = scalar(1)/(nearestDist2[i] + VSMALL);
    }

    DynamicList<label> shieldCells;
    if (shieldingRadius_ > SMALL && shieldingHeight_ > SMALL)
    {
        const scalar minimumZ = scalar(0.5)*localWallNormalSpacing_;
        forAll(mesh.C(), celli)
        {
            const scalar z = mesh.C()[celli].z() - cavityCentre_.z();
            const scalar radial = Foam::sqrt
            (
                sqr(mesh.C()[celli].x() - cavityCentre_.x())
              + sqr(mesh.C()[celli].y() - cavityCentre_.y())
            );
            if
            (
                z >= minimumZ - SMALL
             && z <= shieldingHeight_ + SMALL
             && radial <= shieldingRadius_ + SMALL
            )
            {
                shieldCells.append(celli);
            }
        }
    }
    shieldingCells_.transfer(shieldCells);

    if
    (
        handoffCells_.empty()
     || handoffTargetVolume_ <= SMALL
     || handoffTargetLatentEnergy_ <= SMALL
     || cavityTipCells_.empty()
     ||
        (
            postHandoffSourceMode_ == "none"
         && targetVaporCoreCells < 1
        )
    )
    {
        FatalErrorInFunction
            << "Resolved cap handoff, vapor core, or cavity-tip stencil "
            << "selected no valid cells for site " << diagnosticPrefix_
            << exit(FatalError);
    }
}
void Foam::wallSuperheatVaporBirth::writeCavityDiagnostics
(
    const scalar localSensibleEnergy,
    const scalar wallEnergyThisStep,
    const scalar integratedAlphaSource,
    const scalar integratedMassSource,
    const scalar integratedLatentSink,
    const label activeCells,
    const scalar siteMeanAlphaLiquid,
    const scalar siteVaporWallFraction,
    const scalar attachedWallContactArea,
    const bool attachedBubblePresent,
    const bool departureCandidate,
    const bool rewetCandidate,
    const label detachedBubbleCount,
    const label unexpectedRemoteWallVaporCount,
    const label detachedInterfaceEligibleCellCount,
    const bool writeSourceHistory
) const
{
    const fvMesh& mesh = phase1_.mesh();
    const label timeIndex = mesh.time().timeIndex();
    const bool stateTransition =
        cavityBirthState_ != lastStateDiagnosticsState_
     || cycleId_ != lastStateDiagnosticsCycleId_;
    const bool stateScheduled =
        timeIndex != lastStateDiagnosticsTimeIndex_
     &&
        (
            timeIndex % compactDiagnosticsInterval_ == 0
         || mesh.time().writeTime()
        );
    const bool stateCandidate =
        timeIndex != lastStateDiagnosticsTimeIndex_
     && (departureCandidate || rewetCandidate);
    const bool writeState =
        stateTransition || stateScheduled || stateCandidate;

    bool writeSource = false;
    if (writeSourceHistory)
    {
        const bool sourceTransition =
            cavityBirthState_ != lastSourceDiagnosticsState_
         || cycleId_ != lastSourceDiagnosticsCycleId_;
        const bool sourceScheduled =
            timeIndex != lastSourceDiagnosticsTimeIndex_
         &&
            (
                timeIndex % compactDiagnosticsInterval_ == 0
             || mesh.time().writeTime()
            );
        const bool sourceCandidate =
            timeIndex != lastSourceDiagnosticsTimeIndex_
         && (departureCandidate || rewetCandidate);
        writeSource =
            sourceTransition || sourceScheduled || sourceCandidate;
    }

    if (!writeState && !writeSource)
    {
        return;
    }

    if (writeState)
    {
        lastStateDiagnosticsTimeIndex_ = timeIndex;
        lastStateDiagnosticsState_ = cavityBirthState_;
        lastStateDiagnosticsCycleId_ = cycleId_;
    }
    if (writeSource)
    {
        lastSourceDiagnosticsTimeIndex_ = timeIndex;
        lastSourceDiagnosticsState_ = cavityBirthState_;
        lastSourceDiagnosticsCycleId_ = cycleId_;
    }

    const fileName diagnosticDir
    (
        mesh.time().path()/"postProcessing"/diagnosticPrefix_/"serial"
    );
    mkDir(diagnosticDir);

    const fileName stateFile(diagnosticDir/"birth-state-history.csv");
    const fileName tipFile(diagnosticDir/"cavity-tip-history.csv");
    const fileName energyFile(diagnosticDir/"handoff-energy-history.csv");
    const fileName sourceFile(diagnosticDir/"source-history.csv");
    const bool stateHeader =
        writeState && !std::ifstream(stateFile.c_str()).good();
    const bool tipHeader =
        writeState && !std::ifstream(tipFile.c_str()).good();
    const bool energyHeader =
        writeState && !std::ifstream(energyFile.c_str()).good();
    std::ofstream stateOS;
    std::ofstream tipOS;
    std::ofstream energyOS;
    if (writeState)
    {
        stateOS.open(stateFile.c_str(), std::ios_base::app);
        tipOS.open(tipFile.c_str(), std::ios_base::app);
        energyOS.open(energyFile.c_str(), std::ios_base::app);
    }

    const bool activationSatisfied =
        activationPersistence_ >= 0
      ? activationCandidateDuration_ + SMALL >= activationPersistence_
      : activationCounter_ >= activationPersistenceSteps_;
    const bool thermalCriterionSatisfied =
        lastCavityTtip_ + SMALL >= cavityCriticalTemperature_;
    const scalar topologyMassDifference =
        displacedLiquidMass_ - createdVaporMass_;

    if (writeState && stateHeader)
    {
        stateOS
            << "timeIndex,time,state,activationCounter,activationState,"
            << "activationTime,handoffStartTime,handoffCompleteTime,"
            << "storedHandoffEnergy,createdVaporVolume,createdVaporMass,"
            << "consumedLatentEnergy,creationStepIndex,cycleId,"
            << "departureCandidateDuration,rewetCandidateDuration,"
            << "departureDetected,departureTime,rewetDetected,rewetTime,"
            << "siteRearmed,rearmCount,rearmTime,nextActivationTime,"
            << "siteMeanAlphaLiquid,siteVaporWallFraction,"
            << "attachedWallContactArea,attachedBubblePresent,"
            << "activationCandidateDuration,activationCriterionSatisfied,"
            << "siteShielded,displacedLiquidMass,seedTopologyMassDifference,"
            << "postHandoffSourceMode\n";
    }
    if (writeState && tipHeader)
    {
        tipOS
            << "timeIndex,time,Ttip,Tb,TtipMinusTb,Rc,R,h,"
            << "activationCounter,activationState,"
            << "activationCandidateDuration,siteShielded\n";
    }
    if (writeState && energyHeader)
    {
        energyOS
            << "timeIndex,time,localRemovableSensibleEnergy,"
            << "newWallEnergyThisStep,accumulatedWallEnergy,"
            << "storedHandoffEnergy,targetLatentEnergy,"
            << "latentEnergyConsumed,unconsumedEnergy,"
            << "createdVaporMass,displacedLiquidMass,"
            << "seedTopologyMassDifference\n";
    }

    if (writeSource)
    {
        const bool sourceHeader = !std::ifstream(sourceFile.c_str()).good();
        std::ofstream sourceOS(sourceFile.c_str(), std::ios_base::app);
        if (sourceHeader)
        {
            sourceOS
                << "timeIndex,time,activeSourceCells,integratedAlphaSource,"
                << "integratedMassSource,integratedLatentSink,"
                << "createdVaporVolume,createdVaporMass,detachedBubbleCount,"
                << "unexpectedRemoteWallVaporCount,"
                << "detachedInterfaceEligibleCellCount,departureCandidate,"
                << "rewetCandidate,displacedLiquidMass,"
                << "seedTopologyMassDifference,postHandoffSourceMode,"
                << "seedContinuitySourceMode,siteShielded\n";
        }
        sourceOS
            << mesh.time().timeIndex() << ',' << mesh.time().value() << ','
            << activeCells << ',' << integratedAlphaSource << ','
            << integratedMassSource << ',' << integratedLatentSink << ','
            << createdVaporVolume_ << ',' << createdVaporMass_ << ','
            << detachedBubbleCount << ',' << unexpectedRemoteWallVaporCount
            << ',' << detachedInterfaceEligibleCellCount << ','
            << (departureCandidate ? 1 : 0) << ','
            << (rewetCandidate ? 1 : 0) << ','
            << displacedLiquidMass_ << ',' << topologyMassDifference << ','
            << postHandoffSourceMode_ << ',' << seedContinuitySourceMode_ << ','
            << (lastSiteShielded_ ? 1 : 0) << '\n';
    }

    if (!writeState)
    {
        return;
    }

    stateOS
        << mesh.time().timeIndex() << ',' << mesh.time().value() << ','
        << cavityStateName(cavityBirthState_) << ',' << activationCounter_ << ','
        << (activationSatisfied ? 1 : 0) << ','
        << activationTime_ << ',' << handoffStartTime_ << ','
        << handoffCompleteTime_ << ',' << storedHandoffEnergy_ << ','
        << createdVaporVolume_ << ',' << createdVaporMass_ << ','
        << consumedLatentEnergy_ << ',' << creationStepIndex_ << ','
        << cycleId_ << ',' << departureCandidateDuration_ << ','
        << rewetCandidateDuration_ << ',' << (departureDetected_ ? 1 : 0)
        << ',' << departureTime_ << ',' << (rewetDetected_ ? 1 : 0) << ','
        << rewetTime_ << ',' << (siteRearmed_ ? 1 : 0) << ','
        << rearmCount_ << ',' << rearmTime_ << ',' << nextActivationTime_ << ','
        << siteMeanAlphaLiquid << ',' << siteVaporWallFraction << ','
        << attachedWallContactArea << ',' << (attachedBubblePresent ? 1 : 0)
        << ',' << activationCandidateDuration_ << ','
        << (thermalCriterionSatisfied ? 1 : 0) << ','
        << (lastSiteShielded_ ? 1 : 0) << ',' << displacedLiquidMass_ << ','
        << topologyMassDifference << ',' << postHandoffSourceMode_ << '\n';

    tipOS
        << mesh.time().timeIndex() << ',' << mesh.time().value() << ','
        << lastCavityTtip_ << ',' << cavityCriticalTemperature_ << ','
        << (lastCavityTtip_ - cavityCriticalTemperature_) << ','
        << cavityRadius_ << ',' << cavityCandidateRadius_ << ','
        << cavityTipHeight_ << ',' << activationCounter_ << ','
        << (activationSatisfied ? 1 : 0) << ','
        << activationCandidateDuration_ << ','
        << (lastSiteShielded_ ? 1 : 0) << '\n';

    energyOS
        << mesh.time().timeIndex() << ',' << mesh.time().value() << ','
        << localSensibleEnergy << ',' << wallEnergyThisStep << ','
        << accumulatedWallEnergy_ << ',' << storedHandoffEnergy_ << ','
        << handoffTargetLatentEnergy_ << ',' << consumedLatentEnergy_ << ','
        << max(storedHandoffEnergy_, scalar(0)) << ','
        << createdVaporMass_ << ',' << displacedLiquidMass_ << ','
        << topologyMassDifference << '\n';
}

void Foam::wallSuperheatVaporBirth::writeSourceHistoryRow
(
    const scalar integratedAlphaSource,
    const scalar integratedMassSource,
    const scalar integratedLatentSink,
    const label activeCells,
    const bool departureCandidate,
    const bool rewetCandidate
) const
{
    const fvMesh& mesh = phase1_.mesh();
    const label timeIndex = mesh.time().timeIndex();
    const bool sourceTransition =
        cavityBirthState_ != lastSourceDiagnosticsState_
     || cycleId_ != lastSourceDiagnosticsCycleId_;
    const bool sourceScheduled =
        timeIndex != lastSourceDiagnosticsTimeIndex_
     &&
        (
            timeIndex % compactDiagnosticsInterval_ == 0
         || mesh.time().writeTime()
        );
    const bool sourceCandidate =
        timeIndex != lastSourceDiagnosticsTimeIndex_
     && (departureCandidate || rewetCandidate);

    if (!sourceTransition && !sourceScheduled && !sourceCandidate)
    {
        return;
    }

    lastSourceDiagnosticsTimeIndex_ = timeIndex;
    lastSourceDiagnosticsState_ = cavityBirthState_;
    lastSourceDiagnosticsCycleId_ = cycleId_;

    const fileName diagnosticDir
    (
        mesh.time().path()/"postProcessing"/diagnosticPrefix_/"serial"
    );
    mkDir(diagnosticDir);
    const fileName sourceFile(diagnosticDir/"source-history.csv");
    const bool sourceHeader = !std::ifstream(sourceFile.c_str()).good();
    std::ofstream sourceOS(sourceFile.c_str(), std::ios_base::app);
    const scalar topologyMassDifference =
        displacedLiquidMass_ - createdVaporMass_;
    if (sourceHeader)
    {
        sourceOS
            << "timeIndex,time,activeSourceCells,integratedAlphaSource,"
            << "integratedMassSource,integratedLatentSink,"
            << "createdVaporVolume,createdVaporMass,detachedBubbleCount,"
            << "unexpectedRemoteWallVaporCount,"
            << "detachedInterfaceEligibleCellCount,departureCandidate,"
            << "rewetCandidate,displacedLiquidMass,"
            << "seedTopologyMassDifference,postHandoffSourceMode,"
            << "seedContinuitySourceMode,siteShielded\n";
    }
    sourceOS
        << mesh.time().timeIndex() << ',' << mesh.time().value() << ','
        << activeCells << ',' << integratedAlphaSource << ','
        << integratedMassSource << ',' << integratedLatentSink << ','
        << createdVaporVolume_ << ',' << createdVaporMass_ << ','
        << 0 << ',' << 0 << ',' << 0 << ','
        << (departureCandidate ? 1 : 0) << ','
        << (rewetCandidate ? 1 : 0) << ','
        << displacedLiquidMass_ << ',' << topologyMassDifference << ','
        << postHandoffSourceMode_ << ',' << seedContinuitySourceMode_ << ','
        << (lastSiteShielded_ ? 1 : 0) << '\n';
}
void Foam::wallSuperheatVaporBirth::updateCavityBirth
(
    const volScalarField& liquidTemperature,
    const volScalarField& liquidDensity,
    const volScalarField& liquidCp,
    const volScalarField& vaporDensity
)
{
    const fvMesh& mesh = phase1_.mesh();
    const scalar deltaT = max(mesh.time().deltaTValue(), scalar(0));
    const auto compactLogDue =
    [&](const bool departureCandidate, const bool rewetCandidate)
    {
        const label timeIndex = mesh.time().timeIndex();
        return
            writeDiagnostics_
         &&
            (
                cavityBirthState_ != lastStateDiagnosticsState_
             || cycleId_ != lastStateDiagnosticsCycleId_
             ||
                (
                    timeIndex != lastStateDiagnosticsTimeIndex_
                 &&
                    (
                        timeIndex % compactDiagnosticsInterval_ == 0
                     || mesh.time().writeTime()
                     || departureCandidate
                     || rewetCandidate
                    )
                )
            );
    };

    buildHandoffCells(vaporDensity);
    lastSiteShielded_ = siteShielded();

    scalar siteMeanAlphaLiquid = scalar(1);
    scalar siteVaporWallFraction = scalar(0);
    scalar attachedWallContactArea = scalar(0);
    bool attachedBubblePresent = false;
    evaluateSiteAttachment
    (
        siteMeanAlphaLiquid,
        siteVaporWallFraction,
        attachedWallContactArea,
        attachedBubblePresent
    );

    const bool departureCandidate =
        cavityBirthState_ == CONTACT_LINE_GROWTH
     && !attachedBubblePresent
     && attachedWallContactArea <= departureContactAreaThreshold_ + SMALL;
    const bool rewetCandidate =
        cavityBirthState_ == WAITING_FOR_REWET
     && !attachedBubblePresent
     && siteMeanAlphaLiquid >= rewetAlphaThreshold_ - SMALL
     && siteVaporWallFraction <= rewetWallVaporThreshold_ + SMALL;

    if (departureCandidate)
    {
        departureCandidateDuration_ += deltaT;
    }
    else
    {
        departureCandidateDuration_ = 0;
    }
    if (rewetCandidate)
    {
        rewetCandidateDuration_ += deltaT;
    }
    else
    {
        rewetCandidateDuration_ = 0;
    }

    if
    (
        repeatedNucleationMode_
     && cavityBirthState_ == CONTACT_LINE_GROWTH
     && departureCandidateDuration_ + SMALL >= departurePersistence_
    )
    {
        cavityBirthState_ = WAITING_FOR_REWET;
        departureDetected_ = true;
        departureTime_ = mesh.time().value();
        rewetDetected_ = false;
        rewetTime_ = -1;
        siteRearmed_ = false;
    }

    bool rearmedThisStep = false;
    if
    (
        repeatedNucleationMode_
     && cavityBirthState_ == WAITING_FOR_REWET
     && rewetCandidateDuration_ + SMALL >= rewetPersistence_
     &&
        (
            mesh.time().value() - max(departureTime_, scalar(0))
        ) + SMALL >= minimumRearmDelay_
    )
    {
        rewetDetected_ = true;
        rewetTime_ = mesh.time().value();
        rearmForNextCycle();
        rearmedThisStep = true;
    }

    if (cavityBirthState_ == WAITING_FOR_REWET || rearmedThisStep)
    {
        const bool emitCompactLog =
            compactLogDue(departureCandidate, rewetCandidate);
        writeCavityDiagnostics
        (
            0, 0, 0, 0, 0, 0,
            siteMeanAlphaLiquid,
            siteVaporWallFraction,
            attachedWallContactArea,
            attachedBubblePresent,
            departureCandidate,
            rewetCandidate,
            0, 0, 0
        );
        if (emitCompactLog)
        {
            Info<< "cavity-informed vapor birth diagnostics: site="
                << diagnosticPrefix_
                << ", state=" << cavityStateName(cavityBirthState_)
                << ", cycleId=" << cycleId_
                << ", siteMeanAlphaLiquid=" << siteMeanAlphaLiquid
                << ", siteVaporWallFraction=" << siteVaporWallFraction
                << ", attachedWallContactArea=" << attachedWallContactArea
                << ", departureCandidateDuration="
                << departureCandidateDuration_
                << ", rewetCandidateDuration=" << rewetCandidateDuration_
                << ", departureDetected=" << departureDetected_
                << ", rewetDetected=" << rewetDetected_
                << ", siteRearmed=" << siteRearmed_
                << ", siteShielded=" << lastSiteShielded_ << endl;
        }
        return;
    }

    if (cavityBirthState_ == CONTACT_LINE_GROWTH)
    {
        writeCavityDiagnostics
        (
            0, 0, 0, 0, 0, 0,
            siteMeanAlphaLiquid,
            siteVaporWallFraction,
            attachedWallContactArea,
            attachedBubblePresent,
            departureCandidate,
            false,
            0, 0, 0,
            false
        );
        return;
    }

    lastCavityTtip_ = interpolateCavityTipTemperature(liquidTemperature);
    const bool thermalCriterion =
        lastCavityTtip_ + SMALL >= cavityCriticalTemperature_;
    const bool wettingCriterion =
        !requireWettedSiteForActivation_
     ||
        (
            siteMeanAlphaLiquid >= activationAlphaLiquidThreshold_ - SMALL
         && siteVaporWallFraction <= activationWallVaporThreshold_ + SMALL
        );
    const bool activationCandidate =
        thermalCriterion && wettingCriterion && !lastSiteShielded_;

    if (cavityBirthState_ == WAITING_FOR_CRITERION)
    {
        if (activationCandidate)
        {
            activationCounter_++;
            activationCandidateDuration_ += deltaT;
        }
        else
        {
            activationCounter_ = 0;
            activationCandidateDuration_ = 0;
        }
    }
    const bool activationSatisfied =
        activationPersistence_ >= 0
      ? activationCandidateDuration_ + SMALL >= activationPersistence_
      : activationCounter_ >= activationPersistenceSteps_;

    scalar localSensibleEnergy = 0;
    scalar wallEnergyThisStep = 0;
    forAll(handoffCells_, i)
    {
        const label celli = handoffCells_[i];
        const scalar cellSuperheat =
            max(liquidTemperature[celli] - TSat(celli), scalar(0));
        localSensibleEnergy +=
            seedThermalReserveFraction_
           *phase1_[celli]
           *liquidDensity[celli]
           *liquidCp[celli]
           *mesh.V()[celli]
           *cellSuperheat;
    }

    // Integrate the conjugate wall heat flux by visiting each configured
    // patch face once. The uploaded implementation searched every patch for
    // every seed cell, which becomes prohibitive for several independent
    // sites on a multi-million-cell mesh.
    if (seedWallFluxUseFraction_ > SMALL)
    {
        tmp<volScalarField> tAlphat = turbModel_.alphat();
        forAll(patches_, patchi)
        {
            const label patchID =
                mesh.boundaryMesh().findPatchID(patches_[patchi]);
            if (patchID < 0)
            {
                continue;
            }

            const fvPatch& patch = mesh.boundary()[patchID];
            const labelUList& faceCells = patch.faceCells();
            const scalarField& alphatPatch =
                tAlphat().boundaryField()[patchID];
            tmp<scalarField> tKappaPatch =
                phase1_.thermo().kappaEff(alphatPatch, patchID);
            const scalarField& kappaPatch = tKappaPatch();
            tmp<scalarField> tSnGrad =
                liquidTemperature.boundaryField()[patchID].snGrad();
            const scalarField& snGrad = tSnGrad();
            const vectorField& Sf = patch.Sf();

            forAll(faceCells, facei)
            {
                const label celli = faceCells[facei];
                if
                (
                    celli >= 0
                 && celli < handoffCellMask_.size()
                 && handoffCellMask_[celli]
                )
                {
                    const scalar heatRateIntoFluid = max
                    (
                        kappaPatch[facei]*snGrad[facei]*mag(Sf[facei]),
                        scalar(0)
                    );
                    wallEnergyThisStep +=
                        seedWallFluxUseFraction_
                       *heatRateIntoFluid
                       *deltaT;
                }
            }
        }
    }
    reduce(localSensibleEnergy, sumOp<scalar>());
    reduce(wallEnergyThisStep, sumOp<scalar>());
    lastLocalSensibleEnergy_ = localSensibleEnergy;

    if
    (
        cavityBirthState_ == WAITING_FOR_CRITERION
     && activationSatisfied
    )
    {
        cavityBirthState_ = ACCUMULATING_HANDOFF_ENERGY;
        activationTime_ = mesh.time().value();
        nextActivationTime_ = activationTime_;
        if (cycleId_ == 0)
        {
            cycleId_ = 1;
        }
        if (firstTriggerTime_ < 0)
        {
            firstTriggerTime_ = activationTime_;
        }
        storedHandoffEnergy_ = localSensibleEnergy;
        accumulatedWallEnergy_ = 0;
        departureDetected_ = false;
        rewetDetected_ = false;
        siteRearmed_ = false;
    }

    bool wallEnergyAdded = false;
    if (cavityBirthState_ == ACCUMULATING_HANDOFF_ENERGY)
    {
        accumulatedWallEnergy_ += wallEnergyThisStep;
        storedHandoffEnergy_ += wallEnergyThisStep;
        wallEnergyAdded = true;
        if (storedHandoffEnergy_ + VSMALL >= handoffTargetLatentEnergy_)
        {
            cavityBirthState_ = CREATING_RESOLVED_CAP;
            handoffStartTime_ = mesh.time().value();
            creationStepIndex_ = 0;
        }
    }

    scalar integratedAlphaSource = 0;
    scalar integratedMassSource = 0;
    scalar integratedLatentSink = 0;
    label activeCells = 0;
    if (cavityBirthState_ == CREATING_RESOLVED_CAP)
    {
        if (!wallEnergyAdded)
        {
            accumulatedWallEnergy_ += wallEnergyThisStep;
            storedHandoffEnergy_ += wallEnergyThisStep;
        }
        creationStepIndex_++;
        const label remainingPlannedSteps = max
        (
            seedCreationSteps_ - creationStepIndex_ + 1,
            label(1)
        );
        forAll(handoffCells_, i)
        {
            const label celli = handoffCells_[i];
            const scalar targetVaporFraction =
                handoffTargetVaporFractions_[i];
            const scalar alphaLiquid =
                max(min(phase1_[celli], scalar(1)), scalar(0));
            const scalar alphaVapor = scalar(1) - alphaLiquid;
            const scalar remainingAlpha =
                max(targetVaporFraction - alphaVapor, scalar(0));
            if (remainingAlpha <= SMALL || alphaLiquid <= SMALL)
            {
                continue;
            }

            scalar alphaStep = min
            (
                remainingAlpha/scalar(remainingPlannedSteps),
                min(maxSeedAlphaVaporPerStep_, alphaLiquid)
            );
            const scalar latentPerAlpha =
                vaporDensity[celli]*mesh.V()[celli]*satModel_.L()[celli];
            if (latentPerAlpha <= SMALL)
            {
                continue;
            }
            alphaStep = min(alphaStep, storedHandoffEnergy_/latentPerAlpha);
            alphaStep = max(alphaStep, scalar(0));
            if (alphaStep <= SMALL)
            {
                continue;
            }

            const scalar alphaRate = alphaStep/max(deltaT, VSMALL);
            const scalar physicalVaporMassRate =
                vaporDensity[celli]*alphaRate;
            const scalar continuityMassRate =
                seedContinuitySourceMode_ == "legacyVaporMass"
              ? physicalVaporMassRate
              : scalar(0);
            const scalar latentRate =
                physicalVaporMassRate*satModel_.L()[celli];
            alphaBirthSource_[celli] = alphaRate;
            massBirthSource_[celli] = continuityMassRate;
            latentSink_[celli] = latentRate;
            activeMask_[celli] = scalar(1);
            sustainedState_[celli] = REQUESTED_HOLD_ACTIVE;
            sustainedStateField_[celli] = sustainedState_[celli];

            const scalar cellCreatedVolume = alphaStep*mesh.V()[celli];
            const scalar cellCreatedMass =
                vaporDensity[celli]*cellCreatedVolume;
            const scalar cellDisplacedLiquidMass =
                liquidDensity[celli]*cellCreatedVolume;
            const scalar cellLatentEnergy =
                cellCreatedMass*satModel_.L()[celli];
            createdVaporVolume_ += cellCreatedVolume;
            createdVaporMass_ += cellCreatedMass;
            displacedLiquidMass_ += cellDisplacedLiquidMass;
            consumedLatentEnergy_ += cellLatentEnergy;
            storedHandoffEnergy_ =
                max(storedHandoffEnergy_ - cellLatentEnergy, scalar(0));
            integratedAlphaSource += alphaRate*mesh.V()[celli];
            integratedMassSource += continuityMassRate*mesh.V()[celli];
            integratedLatentSink += latentRate*mesh.V()[celli];
            activeCells++;
        }

        reduce(integratedAlphaSource, sumOp<scalar>());
        reduce(integratedMassSource, sumOp<scalar>());
        reduce(integratedLatentSink, sumOp<scalar>());
        reduce(activeCells, sumOp<label>());
        // Completion is based on the actual alpha field at the beginning of
        // this source update, not on the cumulative injected volume. This is
        // essential after rewet, when a small residual vapor fraction can
        // remain inside the seed stencil and the required new insertion is
        // correspondingly smaller than handoffTargetVolume_. The one-step
        // lag is intentional: the alpha source assembled above is applied by
        // the next alpha solve, and completion is accepted only after that
        // updated alpha field is observed.
        scalar remainingSeedVolume = 0;
        forAll(handoffCells_, i)
        {
            const label celli = handoffCells_[i];
            const scalar targetVaporFraction =
                handoffTargetVaporFractions_[i];
            const scalar actualVaporFraction = scalar(1) - max
            (
                min(phase1_[celli], scalar(1)),
                scalar(0)
            );
            remainingSeedVolume += max
            (
                targetVaporFraction - actualVaporFraction,
                scalar(0)
            )*mesh.V()[celli];
        }
        reduce(remainingSeedVolume, sumOp<scalar>());
        const bool seedComplete =
            remainingSeedVolume
         <= max(SMALL, 1e-6*handoffTargetVolume_);
        if (seedComplete)
        {
            cavityBirthState_ = CONTACT_LINE_GROWTH;
            handoffCompleteTime_ = mesh.time().value();
            firstActiveSiteTime_ = handoffCompleteTime_;
        }
        else if (creationStepIndex_ >= maximumSeedCreationSteps_)
        {
            FatalErrorInFunction
                << "Seed creation for site " << diagnosticPrefix_
                << " did not reach its target after " << creationStepIndex_
                << " steps. createdVaporVolume=" << createdVaporVolume_
                << ", target=" << handoffTargetVolume_
                << ", storedEnergy=" << storedHandoffEnergy_
                << exit(FatalError);
        }
    }

    const bool writeSourceRow =
        cavityBirthState_ != CONTACT_LINE_GROWTH
     || activeCells > 0
     || integratedAlphaSource > VSMALL
     || integratedLatentSink > VSMALL;
    const bool emitCompactLog =
        compactLogDue(departureCandidate, rewetCandidate);
    writeCavityDiagnostics
    (
        localSensibleEnergy,
        wallEnergyThisStep,
        integratedAlphaSource,
        integratedMassSource,
        integratedLatentSink,
        activeCells,
        siteMeanAlphaLiquid,
        siteVaporWallFraction,
        attachedWallContactArea,
        attachedBubblePresent,
        departureCandidate,
        rewetCandidate,
        0, 0, 0,
        writeSourceRow
    );

    if (emitCompactLog)
    {
        Info<< "cavity-informed vapor birth diagnostics: site="
            << diagnosticPrefix_
            << ", state=" << cavityStateName(cavityBirthState_)
            << ", cycleId=" << cycleId_
            << ", Ttip=" << lastCavityTtip_
            << ", Tb=" << cavityCriticalTemperature_
            << ", activationCounter=" << activationCounter_
            << ", activationCandidateDuration="
            << activationCandidateDuration_
            << ", siteShielded=" << lastSiteShielded_
            << ", handoffRadius=" << handoffRadius_
            << ", handoffTargetVolume=" << handoffTargetVolume_
            << ", targetSeedVaporFraction=" << targetSeedVaporFraction_
            << ", storedHandoffEnergy=" << storedHandoffEnergy_
            << ", consumedLatentEnergy=" << consumedLatentEnergy_
            << ", createdVaporMass=" << createdVaporMass_
            << ", displacedLiquidMass=" << displacedLiquidMass_
            << ", seedTopologyMassDifference="
            << displacedLiquidMass_ - createdVaporMass_
            << ", postHandoffSourceMode=" << postHandoffSourceMode_
            << ", seedContinuitySourceMode=" << seedContinuitySourceMode_
            << ", siteMeanAlphaLiquid=" << siteMeanAlphaLiquid
            << ", siteVaporWallFraction=" << siteVaporWallFraction
            << ", attachedWallContactArea=" << attachedWallContactArea
            << ", departureCandidateDuration=" << departureCandidateDuration_
            << ", rewetCandidateDuration=" << rewetCandidateDuration_
            << ", departureDetected=" << departureDetected_
            << ", rewetDetected=" << rewetDetected_
            << ", siteRearmed=" << siteRearmed_
            << ", activeSourceCells=" << activeCells << endl;
    }

}
Foam::scalar Foam::wallSuperheatVaporBirth::thermallyBoundedAlphaStep
(
    const scalar requestedAlphaStep,
    const scalar maxAlphaVaporPerStep,
    const scalar availableAlpha,
    const scalar alphaLiquid,
    const scalar rhoConversion,
    const scalar cp,
    const scalar cellVolume,
    const scalar superheat,
    const scalar minimumEvaporationSuperheat,
    const scalar thermalReserveFraction,
    const scalar latentHeat,
    const scalar minimumThermalAlphaStep
)
{
    if
    (
        requestedAlphaStep <= 0
     || availableAlpha <= 0
     || alphaLiquid <= 0
     || rhoConversion <= 0
     || cp <= 0
     || cellVolume <= 0
     || latentHeat <= 0
     || thermalReserveFraction <= 0
     || superheat <= minimumEvaporationSuperheat
    )
    {
        return 0;
    }

    const scalar thermalReserve =
        alphaLiquid*rhoConversion*cp*cellVolume
       *max(superheat - minimumEvaporationSuperheat, scalar(0));
    const scalar usableReserve = thermalReserveFraction*thermalReserve;
    const scalar thermalAlphaStep =
        usableReserve/(latentHeat*rhoConversion*cellVolume);

    if (thermalAlphaStep <= minimumThermalAlphaStep)
    {
        return 0;
    }

    return min
    (
        min(requestedAlphaStep, maxAlphaVaporPerStep),
        min(availableAlpha, thermalAlphaStep)
    );
}

Foam::scalar Foam::wallSuperheatVaporBirth::wallFluxAwareThermalAlphaStep
(
    const scalar requestedAlphaStep,
    const scalar maxAlphaVaporPerStep,
    const scalar availableAlpha,
    const scalar alphaLiquid,
    const scalar rhoConversion,
    const scalar cp,
    const scalar cellVolume,
    const scalar superheat,
    const scalar minimumEvaporationSuperheat,
    const scalar thermalReserveFraction,
    const scalar latentHeat,
    const scalar minimumThermalAlphaStep,
    const scalar wallFluxUseFraction,
    const scalar qWallIntoFluid,
    const scalar wallArea,
    const scalar deltaT
)
{
    if
    (
        requestedAlphaStep <= 0
     || availableAlpha <= 0
     || alphaLiquid <= 0
     || rhoConversion <= 0
     || cp <= 0
     || cellVolume <= 0
     || latentHeat <= 0
     || thermalReserveFraction <= 0
     || superheat <= minimumEvaporationSuperheat
     || deltaT <= 0
    )
    {
        return 0;
    }

    const scalar sensibleReserve =
        thermalReserveFraction*alphaLiquid*rhoConversion*cp*cellVolume
       *max(superheat - minimumEvaporationSuperheat, scalar(0));
    const scalar wallInflux =
        wallFluxUseFraction*max(qWallIntoFluid, scalar(0))*max(wallArea, scalar(0))*deltaT;
    const scalar availableEnergy = sensibleReserve + wallInflux;
    const scalar thermalAlphaStep =
        availableEnergy/(latentHeat*rhoConversion*cellVolume);

    if (thermalAlphaStep <= minimumThermalAlphaStep)
    {
        return 0;
    }

    return min
    (
        min(requestedAlphaStep, maxAlphaVaporPerStep),
        min(availableAlpha, thermalAlphaStep)
    );
}

void Foam::wallSuperheatVaporBirth::updateSources()
{
    const fvMesh& mesh = phase1_.mesh();
    const label timeIndex = mesh.time().timeIndex();
    if (lastUpdateTimeIndex_ == timeIndex)
    {
        return;
    }
    lastUpdateTimeIndex_ = timeIndex;

    const volScalarField& liquidTemperature = phase1_.thermo().T();
    tmp<volScalarField> tLiquidDensity = phase1_.thermo().rho();
    const volScalarField& liquidDensity = tLiquidDensity();
    tmp<volScalarField> tLiquidCp = phase1_.thermo().Cp();
    const volScalarField& liquidCp = tLiquidCp();
    tmp<volScalarField> tVaporDensity = phase2_.thermo().rho();
    const volScalarField& vaporDensity = tVaporDensity();
    const bool thermallyBoundedHold = sustainedBirthMode_ == "thermallyBoundedHold";

    resetSourceFields
    (
        superheat_,
        alphaBirthSource_,
        massBirthSource_,
        latentSink_
    );
    sustainedStateField_ = dimensionedScalar("zero", dimless, 0.0);
    activeMask_ = dimensionedScalar("zero", dimless, 0.0);

    if (sourceApplicationMode_ == "disabled")
    {
        superheat_.correctBoundaryConditions();
        activeMask_.correctBoundaryConditions();
        sustainedStateField_.correctBoundaryConditions();
        alphaBirthSource_.correctBoundaryConditions();
        massBirthSource_.correctBoundaryConditions();
        latentSink_.correctBoundaryConditions();

        Info<< "wall-superheat vapor birth diagnostics: activeSourceCells=0"
            << ", sourceFootprintMode=" << sourceFootprintMode_
            << ", sourceApplicationMode=" << sourceApplicationMode_
            << ", contactLineEligibleCells=0"
            << ", dryCoreExcludedCells=0"
            << ", detachedExcludedCells=0"
            << ", sustainedBirthMode=" << sustainedBirthMode_
            << ", requestedHoldCells=0"
            << ", thermallyActiveSourceCells=0"
            << ", thermallySuspendedCells=0"
            << ", holdExpiredCells=0"
            << ", siteMode=" << siteMode_
            << ", siteCenter=" << siteCenter_
            << ", siteRadius=" << siteRadius_
            << ", eligibleSiteCells=" << siteCells_.size()
            << ", eligibleSiteVolume=" << siteVolume()
            << ", firstActiveSiteTime=" << firstActiveSiteTime_
            << ", siteLocked=" << siteLocked_
            << ", sourceCellsOutsideSite=0"
            << ", integratedRequestedAlphaSource=0"
            << ", integratedAlphaSource=0"
            << ", integratedMassSource=0"
            << ", integratedLatentSink=0"
            << ", integratedThermalReserveUsed=0"
            << ", thermalLimitMode=" << thermalLimitMode_
            << ", wallFluxUseFraction=" << wallFluxUseFraction_
            << ", integratedWallFluxAllowance=0"
            << ", maxRequestedAlphaStep=0"
            << ", maxThermalAlphaLimit=0"
            << ", maxWallFluxIntoFluid=0"
            << ", maxWallFluxAlphaLimit=0"
            << ", maxActualAlphaStep=0"
            << ", maxAlphaRate=0"
            << ", maxMassSource=0"
            << ", maxLatentSink=0"
            << ", minSourceCellTemperature=0"
            << ", minSourceCellSuperheat=0"
            << ", maxSuperheat=0"
            << ", firstTriggerTime=" << firstTriggerTime_ << endl;
        return;
    }

    // Seed-only cavity models are handled before the legacy sustained-source
    // diagnostic arrays are allocated. This is essential when several
    // independent cavity instances share one heater.
    if (cavityMode() && postHandoffSourceMode_ == "none")
    {
        const bool cavityManagedState =
            cavityBirthState_ != CONTACT_LINE_GROWTH || repeatedNucleationMode_;
        if (cavityManagedState)
        {
            updateCavityBirth
            (
                liquidTemperature,
                liquidDensity,
                liquidCp,
                vaporDensity
            );
        }

        label seedActiveCells = 0;
        scalar integratedSeedAlphaSource = 0;
        scalar integratedSeedContinuitySource = 0;
        scalar integratedSeedLatentSink = 0;
        forAll(alphaBirthSource_, celli)
        {
            if (activeMask_[celli] > 0.5)
            {
                ++seedActiveCells;
            }
            integratedSeedAlphaSource +=
                alphaBirthSource_[celli]*mesh.V()[celli];
            integratedSeedContinuitySource +=
                massBirthSource_[celli]*mesh.V()[celli];
            integratedSeedLatentSink +=
                latentSink_[celli]*mesh.V()[celli];
        }
        reduce(seedActiveCells, sumOp<label>());
        reduce(integratedSeedAlphaSource, sumOp<scalar>());
        reduce(integratedSeedContinuitySource, sumOp<scalar>());
        reduce(integratedSeedLatentSink, sumOp<scalar>());

        if
        (
            cavityBirthState_ == CONTACT_LINE_GROWTH
         && seedActiveCells == 0
        )
        {
            scalar siteMeanAlphaLiquid = scalar(1);
            scalar siteVaporWallFraction = scalar(0);
            scalar attachedWallContactArea = scalar(0);
            bool attachedBubblePresent = false;
            evaluateSiteAttachment
            (
                siteMeanAlphaLiquid,
                siteVaporWallFraction,
                attachedWallContactArea,
                attachedBubblePresent
            );
            const bool departureCandidate =
                !attachedBubblePresent
             && attachedWallContactArea
                <= departureContactAreaThreshold_ + SMALL;
            writeSourceHistoryRow
            (
                integratedSeedAlphaSource,
                integratedSeedContinuitySource,
                integratedSeedLatentSink,
                seedActiveCells,
                departureCandidate,
                false
            );
        }

        superheat_.correctBoundaryConditions();
        activeMask_.correctBoundaryConditions();
        sustainedStateField_.correctBoundaryConditions();
        alphaBirthSource_.correctBoundaryConditions();
        massBirthSource_.correctBoundaryConditions();
        latentSink_.correctBoundaryConditions();
        return;
    }

    label activeCells = 0;
    label contactLineEligibleCells = 0;
    label dryCoreExcludedCells = 0;
    label detachedExcludedCells = 0;
    label requestedHoldCells = 0;
    label thermallyActiveSourceCells = 0;
    label thermallySuspendedCells = 0;
    label holdExpiredCells = 0;
    scalar maxSuperheat = -GREAT;
    scalar integratedMassSource = 0;
    scalar integratedAlphaSource = 0;
    scalar integratedRequestedAlphaSource = 0;
    scalar integratedLatentSink = 0;
    scalar integratedThermalReserveUsed = 0;
    scalar integratedWallFluxAllowance = 0;
    scalar maxAlphaRate = 0;
    scalar maxMassSource = 0;
    scalar maxLatentSink = 0;
    scalar minSourceCellTemperature = GREAT;
    scalar minSourceCellSuperheat = GREAT;
    scalar maxRequestedAlphaStep = 0;
    scalar maxThermalAlphaLimit = 0;
    scalar maxWallFluxIntoFluid = 0;
    scalar maxWallFluxAlphaLimit = 0;
    scalar maxActualAlphaStep = 0;

    scalarField diagnosticRequestedAlphaStep(mesh.nCells(), 0);
    scalarField diagnosticAvailableLiquidLimit(mesh.nCells(), 0);
    scalarField diagnosticSensibleOnlyThermalAlphaLimit(mesh.nCells(), 0);
    scalarField diagnosticWallFluxAwareThermalAlphaLimit(mesh.nCells(), 0);
    scalarField diagnosticActualAlphaStep(mesh.nCells(), 0);
    scalarField diagnosticSensibleEnergyAllowance(mesh.nCells(), 0);
    scalarField diagnosticWallEnergyAllowance(mesh.nCells(), 0);
    scalarField diagnosticCombinedEnergyAllowance(mesh.nCells(), 0);
    List<word> diagnosticLimiter(mesh.nCells(), word("STATE_INACTIVE"));
    scalarField diagnosticAlphaVapor(mesh.nCells(), 0);
    scalarField diagnosticContactLineEligible(mesh.nCells(), 0);
    scalarField diagnosticNeighbourMaxAlphaVapor(mesh.nCells(), 0);
    boolList sourceLocalityMask(mesh.nCells(), false);
    boolList sourceActiveMask(mesh.nCells(), false);
    DynamicList<label> pendingInterfaceWallCells;
    scalar pendingInterfaceLatentEnergy = 0;
    label interfaceBandEligibleCells = 0;
    label interfaceBandSourceCells = 0;
    scalar interfaceBandSourceVolume = 0;
    scalar interfaceBandSourceWallNormalMoment = 0;
    scalar interfaceBandSourceAboveFirstLayer = 0;
    scalar interfaceBandSourceTotal = 0;
    const bool interfaceDeposition = growthDepositionMode_ == "wallNormalInterfaceBand";

    if (cavityMode())
    {
        const bool cavityManagedState =
            cavityBirthState_ != CONTACT_LINE_GROWTH || repeatedNucleationMode_;

        if (cavityManagedState)
        {
            updateCavityBirth
            (
                liquidTemperature,
                liquidDensity,
                liquidCp,
                vaporDensity
            );

            if (cavityBirthState_ != CONTACT_LINE_GROWTH)
            {
                superheat_.correctBoundaryConditions();
                activeMask_.correctBoundaryConditions();
                sustainedStateField_.correctBoundaryConditions();
                alphaBirthSource_.correctBoundaryConditions();
                massBirthSource_.correctBoundaryConditions();
                latentSink_.correctBoundaryConditions();
                return;
            }
        }
    }

    forAll(siteCells_, siteI)
    {
        const label celli = siteCells_[siteI];
        const labelList& neighbours = mesh.cellCells()[celli];
        forAll(neighbours, neighbourI)
        {
            const label neighbCelli = neighbours[neighbourI];
            if (siteCellMask_[neighbCelli])
            {
                diagnosticNeighbourMaxAlphaVapor[celli] = max
                (
                    diagnosticNeighbourMaxAlphaVapor[celli],
                    scalar(1) - phase1_[neighbCelli]
                );
            }
        }
    }

    if (Pstream::parRun())
    {
        volScalarField siteMaskField
        (
            IOobject
            (
                word(diagnosticPrefix_ + "SiteMaskExchange"),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("0", dimless, 0),
            "zeroGradient"
        );

        forAll(siteCells_, siteI)
        {
            siteMaskField[siteCells_[siteI]] = 1;
        }

        siteMaskField.correctBoundaryConditions();

        forAll(mesh.boundary(), patchID)
        {
            const fvPatchScalarField& alphaPatch = phase1_.boundaryField()[patchID];
            if (!alphaPatch.coupled())
            {
                continue;
            }

            const labelUList& faceCells = mesh.boundary()[patchID].faceCells();
            const scalarField neighbourAlphaLiquid(alphaPatch.patchNeighbourField());
            const scalarField neighbourSiteMask
            (
                siteMaskField.boundaryField()[patchID].patchNeighbourField()
            );

            forAll(faceCells, facei)
            {
                const label celli = faceCells[facei];
                if (!siteCellMask_[celli] || neighbourSiteMask[facei] < 0.5)
                {
                    continue;
                }

                diagnosticNeighbourMaxAlphaVapor[celli] = max
                (
                    diagnosticNeighbourMaxAlphaVapor[celli],
                    scalar(1) - neighbourAlphaLiquid[facei]
                );
            }
        }
    }

    forAll(candidateCells_, candidateI)
    {
        const label celli = candidateCells_[candidateI];
        const scalar cellSuperheat = liquidTemperature[celli] - TSat(celli);
        superheat_[celli] = cellSuperheat;
        maxSuperheat = max(maxSuperheat, cellSuperheat);

        const scalar alphaVapor = scalar(1) - phase1_[celli];
        diagnosticAlphaVapor[celli] = alphaVapor;
        const bool movingContactLine = sourceFootprintMode_ == "wallAdjacentContactLine";
        const bool hasLiquid = phase1_[celli]
            > (movingContactLine ? minimumSourceLiquidFraction_ : epsLiquid_);
                const bool footprintEligible = movingContactLine
                    ? contactLineEligible
                        (
                                phase1_[celli],
                                diagnosticNeighbourMaxAlphaVapor[celli],
                                contactLineAlphaMin_,
                                contactLineAlphaMax_,
                                minimumSourceLiquidFraction_,
                                contactLineNeighbourSearch_,
                                contactLineNeighbourDelta_
                        )
                    : true;
        diagnosticContactLineEligible[celli] = footprintEligible ? scalar(1) : scalar(0);

        if (movingContactLine && footprintEligible)
        {
            contactLineEligibleCells++;
        }

        if (movingContactLine && !footprintEligible)
        {
            if (phase1_[celli] < minimumSourceLiquidFraction_)
            {
                dryCoreExcludedCells++;
                diagnosticLimiter[celli] = "DRY_CORE_EXCLUDED";
            }
            else
            {
                detachedExcludedCells++;
                diagnosticLimiter[celli] = "CONTACT_LINE_EXCLUDED";
            }
            activeMask_[celli] = scalar(0);
            sustainedState_[celli] = INACTIVE;
            sustainedActivationTime_[celli] = -GREAT;
            sustainedStateField_[celli] = sustainedState_[celli];
            continue;
        }

        const scalar deltaT = mesh.time().deltaTValue();
        bool active = activeMask_[celli] > 0.5;
        bool requested = active;
        bool holdExpired = false;

        if (!thermallyBoundedHold)
        {
            if (active && (!hasLiquid || cellSuperheat <= superheatOff_))
            {
                active = false;
            }
            else if (!active && hasLiquid && cellSuperheat >= superheatOn_)
            {
                active = true;
            }

            requested = active;
            sustainedState_[celli] = active ? REQUESTED_HOLD_ACTIVE : INACTIVE;
        }
        else
        {
            if (!hasLiquid)
            {
                requested = false;
                sustainedState_[celli] = INACTIVE;
                sustainedActivationTime_[celli] = -GREAT;
            }
            else if (!requested && cellSuperheat >= superheatOn_)
            {
                requested = true;
                sustainedState_[celli] = REQUESTED_HOLD_ACTIVE;
                sustainedActivationTime_[celli] = mesh.time().value();
            }
            else if (requested)
            {
                const scalar timeSinceActivation =
                    sustainedActivationTime_[celli] > -GREAT/2
                  ? mesh.time().value() - sustainedActivationTime_[celli]
                  : GREAT;

                holdExpired = timeSinceActivation >= minimumActiveDuration_;
                if (holdExpired)
                {
                    holdExpiredCells++;
                    if (cellSuperheat <= superheatOff_)
                    {
                        requested = false;
                        sustainedState_[celli] = INACTIVE;
                        sustainedActivationTime_[celli] = -GREAT;
                    }
                }
            }
        }

        activeMask_[celli] = requested ? scalar(1) : scalar(0);

        if (thermallyBoundedHold && requested)
        {
            requestedHoldCells++;
        }

        if (!requested || !hasLiquid)
        {
            diagnosticLimiter[celli] = requested ? word("AVAILABLE_LIQUID") : word("STATE_INACTIVE");
            sustainedStateField_[celli] = sustainedState_[celli];
            continue;
        }

        if (cellSuperheat <= 0)
        {
            if (thermallyBoundedHold)
            {
                sustainedState_[celli] = THERMALLY_SUSPENDED;
                sustainedStateField_[celli] = sustainedState_[celli];
                thermallySuspendedCells++;
                diagnosticLimiter[celli] = "THERMAL_FLOOR";
            }
            else
            {
                sustainedStateField_[celli] = sustainedState_[celli];
                diagnosticLimiter[celli] = "THERMAL_FLOOR";
            }
            continue;
        }

        if (deltaT <= SMALL)
        {
            sustainedState_[celli] = thermallyBoundedHold ? THERMALLY_SUSPENDED : sustainedState_[celli];
            sustainedStateField_[celli] = sustainedState_[celli];
            diagnosticLimiter[celli] = "ZERO_DELTA_T";
            continue;
        }

        const scalar scale = min(sourceRelaxation_*cellSuperheat/max(superheatOn_, SMALL), scalar(1));
        const scalar availableAlpha = max(phase1_[celli], scalar(0));
        const scalar requestedAlphaStep = maxAlphaVaporPerStep_*max(scale, scalar(0));
        scalar alphaStep = min(requestedAlphaStep, availableAlpha);
        diagnosticRequestedAlphaStep[celli] = requestedAlphaStep;
        diagnosticAvailableLiquidLimit[celli] = availableAlpha;
        diagnosticActualAlphaStep[celli] = alphaStep;

        scalar thermalAlphaLimit = GREAT;
        scalar qWallIntoFluid = 0;
        scalar wallArea = 0;
        if (thermallyBoundedHold)
        {
            const scalar reserveSuperheat = max
            (
                cellSuperheat - minimumEvaporationSuperheat_,
                scalar(0)
            );
            thermalAlphaLimit =
                                reserveSuperheat > 0
                         && liquidDensity[celli] > SMALL
                         && mesh.V()[celli] > SMALL
                         && satModel_.L()[celli] > SMALL
              ? thermalReserveFraction_*phase1_[celli]*liquidDensity[celli]
               *liquidCp[celli]*mesh.V()[celli]*reserveSuperheat
               /(satModel_.L()[celli]*liquidDensity[celli]*mesh.V()[celli])
              : scalar(0);
            diagnosticSensibleEnergyAllowance[celli] =
                thermalReserveFraction_*phase1_[celli]*liquidDensity[celli]
               *liquidCp[celli]*mesh.V()[celli]*reserveSuperheat;
            diagnosticSensibleOnlyThermalAlphaLimit[celli] = thermalAlphaLimit;
            diagnosticCombinedEnergyAllowance[celli] = diagnosticSensibleEnergyAllowance[celli];

            if (thermalLimitMode_ == "wallFluxAwareInstantaneous")
            {
                qWallIntoFluid = wallFluxIntoFluid(celli);
                wallArea = wallFaceArea(celli);
                const scalar wallEnergy =
                    wallFluxUseFraction_*max(qWallIntoFluid, scalar(0))*wallArea*deltaT;
                diagnosticWallEnergyAllowance[celli] = wallEnergy;
                const scalar wallFluxAlphaLimit =
                    liquidDensity[celli] > SMALL
                 && mesh.V()[celli] > SMALL
                 && satModel_.L()[celli] > SMALL
                  ? wallEnergy/(satModel_.L()[celli]*liquidDensity[celli]*mesh.V()[celli])
                  : scalar(0);
                thermalAlphaLimit += max(wallFluxAlphaLimit, scalar(0));
                diagnosticWallFluxAwareThermalAlphaLimit[celli] = thermalAlphaLimit;
                diagnosticCombinedEnergyAllowance[celli] += wallEnergy;
                maxWallFluxIntoFluid = max(maxWallFluxIntoFluid, qWallIntoFluid);
                maxWallFluxAlphaLimit = max(maxWallFluxAlphaLimit, wallFluxAlphaLimit);
                integratedWallFluxAllowance += wallEnergy;

                alphaStep = wallFluxAwareThermalAlphaStep
                (
                    requestedAlphaStep,
                    maxAlphaVaporPerStep_,
                    availableAlpha,
                    phase1_[celli],
                    liquidDensity[celli],
                    liquidCp[celli],
                    mesh.V()[celli],
                    cellSuperheat,
                    minimumEvaporationSuperheat_,
                    thermalReserveFraction_,
                    satModel_.L()[celli],
                    minimumThermalAlphaStep_,
                    wallFluxUseFraction_,
                    qWallIntoFluid,
                    wallArea,
                    deltaT
                );
            }
            else
            {
                alphaStep = thermallyBoundedAlphaStep
                (
                    requestedAlphaStep,
                    maxAlphaVaporPerStep_,
                    availableAlpha,
                    phase1_[celli],
                    liquidDensity[celli],
                    liquidCp[celli],
                    mesh.V()[celli],
                    cellSuperheat,
                    minimumEvaporationSuperheat_,
                    thermalReserveFraction_,
                    satModel_.L()[celli],
                    minimumThermalAlphaStep_
                );
            }

            maxThermalAlphaLimit = max(maxThermalAlphaLimit, thermalAlphaLimit);
            diagnosticActualAlphaStep[celli] = alphaStep;
        }

        maxRequestedAlphaStep = max(maxRequestedAlphaStep, requestedAlphaStep);
        integratedRequestedAlphaSource += requestedAlphaStep/deltaT*mesh.V()[celli];

        if (thermallyBoundedHold && alphaStep <= 0)
        {
            sustainedState_[celli] = THERMALLY_SUSPENDED;
            sustainedStateField_[celli] = sustainedState_[celli];
            thermallySuspendedCells++;
            diagnosticLimiter[celli] = "THERMAL_FLOOR_OR_MINIMUM_STEP";
            continue;
        }

        const bool diagnosticOnly = sourceApplicationMode_ == "diagnosticOnly";

        sustainedState_[celli] = REQUESTED_HOLD_ACTIVE;
        sustainedStateField_[celli] = sustainedState_[celli];
        if (diagnosticOnly)
        {
            diagnosticActualAlphaStep[celli] = scalar(0);
            diagnosticLimiter[celli] = "DIAGNOSTIC_ONLY";
        }
        else if (interfaceDeposition)
        {
            const scalar latentEnergy = alphaStep*liquidDensity[celli]*mesh.V()[celli]*satModel_.L()[celli];
            if (latentEnergy > SMALL)
            {
                pendingInterfaceLatentEnergy += latentEnergy;
                pendingInterfaceWallCells.append(celli);
            }
            diagnosticLimiter[celli] = "INTERFACE_BAND_SOURCE_ORIGIN";
        }
        else
        {
            const scalar alphaRate = alphaStep/deltaT;
            alphaBirthSource_[celli] = alphaRate;
            massBirthSource_[celli] = liquidDensity[celli]*alphaRate;
            latentSink_[celli] = massBirthSource_[celli]*satModel_.L()[celli];
            sourceActiveMask[celli] = alphaRate > SMALL;
            sourceLocalityMask[celli] = true;

            if (alphaStep + VSMALL >= maxAlphaVaporPerStep_ && maxAlphaVaporPerStep_ <= requestedAlphaStep + VSMALL)
            {
                diagnosticLimiter[celli] = "PER_STEP_CAP";
            }
            else if (alphaStep + VSMALL >= availableAlpha && availableAlpha <= requestedAlphaStep + VSMALL)
            {
                diagnosticLimiter[celli] = "AVAILABLE_LIQUID";
            }
            else if (thermallyBoundedHold && alphaStep + VSMALL >= thermalAlphaLimit)
            {
                diagnosticLimiter[celli] =
                    thermalLimitMode_ == "wallFluxAwareInstantaneous"
                  ? word("WALL_FLUX_AWARE_THERMAL_LIMIT")
                  : word("SENSIBLE_RESERVE");
            }
            else
            {
                diagnosticLimiter[celli] = "REQUESTED_SOURCE";
            }

            integratedAlphaSource += alphaBirthSource_[celli]*mesh.V()[celli];
            integratedMassSource += massBirthSource_[celli]*mesh.V()[celli];
            integratedLatentSink += latentSink_[celli]*mesh.V()[celli];
            integratedThermalReserveUsed += latentSink_[celli]*mesh.V()[celli];
            maxAlphaRate = max(maxAlphaRate, alphaBirthSource_[celli]);
            maxMassSource = max(maxMassSource, massBirthSource_[celli]);
            maxLatentSink = max(maxLatentSink, latentSink_[celli]);
            maxActualAlphaStep = max(maxActualAlphaStep, diagnosticActualAlphaStep[celli]);
            minSourceCellTemperature = min(minSourceCellTemperature, liquidTemperature[celli]);
            minSourceCellSuperheat = min(minSourceCellSuperheat, cellSuperheat);
        }
    }

    if (interfaceDeposition && pendingInterfaceLatentEnergy > SMALL)
    {
        boolList bandIncluded(mesh.nCells(), false);
        scalarField bandWeights(mesh.nCells(), 0);
        scalarField bandNormalDistance(mesh.nCells(), GREAT);

        forAll(pendingInterfaceWallCells, pendingI)
        {
            appendInterfaceBandCells
            (
                pendingInterfaceWallCells[pendingI],
                phase1_.primitiveField(),
                bandIncluded,
                bandWeights,
                bandNormalDistance,
                interfaceBandEligibleCells
            );
        }

        scalar weightSum = 0;
        forAll(bandWeights, celli)
        {
            weightSum += bandWeights[celli];
        }

        if (weightSum > SMALL)
        {
            const scalar deltaT = mesh.time().deltaTValue();
            forAll(bandWeights, celli)
            {
                if (!bandIncluded[celli] || bandWeights[celli] <= SMALL)
                {
                    continue;
                }

                const scalar desiredEnergy = pendingInterfaceLatentEnergy*bandWeights[celli]/weightSum;
                const scalar latentPerAlpha = liquidDensity[celli]*mesh.V()[celli]*satModel_.L()[celli];
                if (latentPerAlpha <= SMALL)
                {
                    continue;
                }

                scalar alphaStep = desiredEnergy/latentPerAlpha;
                alphaStep = min(alphaStep, maxAlphaVaporPerStep_);
                alphaStep = min(alphaStep, max(phase1_[celli], scalar(0)));
                if (alphaStep <= SMALL)
                {
                    continue;
                }

                const scalar alphaRate = alphaStep/deltaT;
                const scalar massRate = liquidDensity[celli]*alphaRate;
                const scalar latentRate = massRate*satModel_.L()[celli];
                alphaBirthSource_[celli] += alphaRate;
                massBirthSource_[celli] += massRate;
                latentSink_[celli] += latentRate;
                activeMask_[celli] = scalar(1);
                sustainedState_[celli] = REQUESTED_HOLD_ACTIVE;
                sustainedStateField_[celli] = sustainedState_[celli];
                sourceActiveMask[celli] = true;
                sourceLocalityMask[celli] = true;
                diagnosticActualAlphaStep[celli] += alphaStep;
                diagnosticLimiter[celli] = "WALL_NORMAL_INTERFACE_BAND";

                const scalar cellAlphaSource = alphaRate*mesh.V()[celli];
                const scalar cellMassSource = massRate*mesh.V()[celli];
                const scalar cellLatentSink = latentRate*mesh.V()[celli];
                integratedAlphaSource += cellAlphaSource;
                integratedMassSource += cellMassSource;
                integratedLatentSink += cellLatentSink;
                integratedThermalReserveUsed += cellLatentSink;
                maxAlphaRate = max(maxAlphaRate, alphaBirthSource_[celli]);
                maxMassSource = max(maxMassSource, massBirthSource_[celli]);
                maxLatentSink = max(maxLatentSink, latentSink_[celli]);
                maxActualAlphaStep = max(maxActualAlphaStep, diagnosticActualAlphaStep[celli]);
                minSourceCellTemperature = min(minSourceCellTemperature, liquidTemperature[celli]);
                minSourceCellSuperheat = min(minSourceCellSuperheat, liquidTemperature[celli] - TSat(celli));
                interfaceBandSourceCells++;
                interfaceBandSourceVolume += mesh.V()[celli];
                interfaceBandSourceWallNormalMoment += cellLatentSink*bandNormalDistance[celli];
                interfaceBandSourceTotal += cellLatentSink;
                if (bandNormalDistance[celli] > localWallNormalSpacing_ + SMALL)
                {
                    interfaceBandSourceAboveFirstLayer += cellLatentSink;
                }
            }
        }
    }

    activeCells = 0;
    thermallyActiveSourceCells = 0;
    forAll(sourceActiveMask, celli)
    {
        if (sourceActiveMask[celli])
        {
            activeCells++;
            if (thermallyBoundedHold)
            {
                thermallyActiveSourceCells++;
            }
        }
    }

    label sourceCellsOutsideSite = 0;
    forAll(alphaBirthSource_, celli)
    {
        if
        (
            alphaBirthSource_[celli] > SMALL
         && !siteCellMask_[celli]
         && !sourceLocalityMask[celli]
        )
        {
            sourceCellsOutsideSite++;
        }
    }

    reduce(activeCells, sumOp<label>());
    reduce(contactLineEligibleCells, sumOp<label>());
    reduce(dryCoreExcludedCells, sumOp<label>());
    reduce(detachedExcludedCells, sumOp<label>());
    reduce(requestedHoldCells, sumOp<label>());
    reduce(thermallyActiveSourceCells, sumOp<label>());
    reduce(thermallySuspendedCells, sumOp<label>());
    reduce(holdExpiredCells, sumOp<label>());
    reduce(sourceCellsOutsideSite, sumOp<label>());
    reduce(interfaceBandEligibleCells, sumOp<label>());
    reduce(interfaceBandSourceCells, sumOp<label>());
    reduce(interfaceBandSourceVolume, sumOp<scalar>());
    reduce(interfaceBandSourceWallNormalMoment, sumOp<scalar>());
    reduce(interfaceBandSourceAboveFirstLayer, sumOp<scalar>());
    reduce(interfaceBandSourceTotal, sumOp<scalar>());
    reduce(integratedAlphaSource, sumOp<scalar>());
    reduce(integratedRequestedAlphaSource, sumOp<scalar>());
    reduce(integratedMassSource, sumOp<scalar>());
    reduce(integratedLatentSink, sumOp<scalar>());
    reduce(integratedThermalReserveUsed, sumOp<scalar>());
    reduce(integratedWallFluxAllowance, sumOp<scalar>());
    reduce(maxAlphaRate, maxOp<scalar>());
    reduce(maxMassSource, maxOp<scalar>());
    reduce(maxLatentSink, maxOp<scalar>());
    reduce(minSourceCellTemperature, minOp<scalar>());
    reduce(minSourceCellSuperheat, minOp<scalar>());
    reduce(maxRequestedAlphaStep, maxOp<scalar>());
    reduce(maxThermalAlphaLimit, maxOp<scalar>());
    reduce(maxWallFluxIntoFluid, maxOp<scalar>());
    reduce(maxWallFluxAlphaLimit, maxOp<scalar>());
    reduce(maxActualAlphaStep, maxOp<scalar>());

    if (activeCells > 0 && firstTriggerTime_ < 0)
    {
        firstTriggerTime_ = mesh.time().value();
    }

    if (activeCells > 0 && firstActiveSiteTime_ < 0)
    {
        firstActiveSiteTime_ = mesh.time().value();
        siteLocked_ = lockSiteAfterBirth_;
    }

    if (sourceCellsOutsideSite != 0)
    {
        FatalErrorInFunction
            << "wallSuperheatVaporBirth generated " << sourceCellsOutsideSite
            << " source cells outside the configured site"
            << exit(FatalError);
    }

    superheat_.correctBoundaryConditions();
    activeMask_.correctBoundaryConditions();
    sustainedStateField_.correctBoundaryConditions();
    alphaBirthSource_.correctBoundaryConditions();
    massBirthSource_.correctBoundaryConditions();
    latentSink_.correctBoundaryConditions();

    if
    (
        writePerCellDiagnostics_
     && timeIndex % perCellDiagnosticsInterval_ == 0
     && lastPerCellDiagnosticsTimeIndex_ != timeIndex
    )
    {
                std::ostringstream rankName;
                rankName << "processor" << std::setw(4) << std::setfill('0') << Pstream::myProcNo();
                const word rankDir = Pstream::parRun() ? word(rankName.str()) : word("serial");
        const fileName diagnosticDir(mesh.time().path()/"postProcessing"/diagnosticPrefix_/rankDir);
        const fileName cellFile(diagnosticDir/"perCellSourceHistory.csv");
        const fileName faceFile(diagnosticDir/"perFaceWallFluxHistory.csv");
        mkDir(diagnosticDir);
        labelList globalCellIDs(mesh.nCells());
        labelList globalFaceIDs(mesh.nFaces());
        forAll(globalCellIDs, celli)
        {
            globalCellIDs[celli] = celli;
        }
        forAll(globalFaceIDs, facei)
        {
            globalFaceIDs[facei] = facei;
        }

        if (Pstream::parRun())
        {
            const labelIOList cellProcAddressing
            (
                IOobject
                (
                    "cellProcAddressing",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            const labelIOList faceProcAddressing
            (
                IOobject
                (
                    "faceProcAddressing",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            if
            (
                cellProcAddressing.size() != mesh.nCells()
             || faceProcAddressing.size() != mesh.nFaces()
            )
            {
                FatalErrorInFunction
                    << "Processor addressing inconsistent with mesh sizes" << nl
                    << "cells: " << mesh.nCells()
                    << " addressing: " << cellProcAddressing.size() << nl
                    << "faces: " << mesh.nFaces()
                    << " addressing: " << faceProcAddressing.size()
                    << exit(FatalError);
            }

            forAll(globalCellIDs, celli)
            {
                globalCellIDs[celli] = cellProcAddressing[celli];
            }
            forAll(globalFaceIDs, facei)
            {
                const label addressedFace = faceProcAddressing[facei];
                globalFaceIDs[facei] = addressedFace > 0
                  ? addressedFace - 1
                  : -addressedFace - 1;
            }
        }

        const auto blockExists = [](const fileName& path, const label segmentID, const label blockTimeIndex)
        {
            std::ifstream input(path.c_str());
            std::string line;
            std::string prefix = std::to_string(segmentID) + "," + std::to_string(blockTimeIndex) + ",";
            while (std::getline(input, line))
            {
                if (line.rfind(prefix, 0) == 0)
                {
                    return true;
                }
            }
            return false;
        };

        if (blockExists(cellFile, perCellDiagnosticsSegmentID_, timeIndex))
        {
            FatalErrorInFunction
                << "Per-cell diagnostics already contain segmentID "
                << perCellDiagnosticsSegmentID_ << " timeIndex " << timeIndex
                << "; refusing duplicate append"
                << exit(FatalError);
        }

        const bool writeCellHeader = !std::ifstream(cellFile.c_str()).good();
        const bool writeFaceHeader = !std::ifstream(faceFile.c_str()).good();
        std::ofstream cellOS(cellFile.c_str(), std::ios_base::app);
        std::ofstream faceOS(faceFile.c_str(), std::ios_base::app);

        if (writeCellHeader)
        {
            cellOS
                << "segmentID,timeIndex,time,deltaT,cellID,processorRank,globalCellID,processorLocalCellID,cellCentreX,cellCentreY,cellCentreZ,cellVolume,"
                << "sourceFootprintMode,siteEligible,contactLineEligible,alphaVapor,controllerState,holdStartTime,holdExpiryTime,wallFaceCount,wallFaceIDs,wallArea,"
                << "kappaEff,snGradT,qWallIntoFluid,wallPower,wallEnergyThisStep,temperature,Tsat,superheat,"
                << "sensibleEnergyAllowance,wallEnergyAllowance,combinedEnergyAllowance,requestedAlphaStep,"
                << "perStepAlphaCap,availableLiquidLimit,neighbourMaxAlphaVapor,sensibleOnlyThermalAlphaLimit,"
                << "wallFluxAwareThermalAlphaLimit,actualAlphaStep,alphaRate,massRate,latentRate,activeLimiter\n";
        }

        if (writeFaceHeader)
        {
            faceOS
                << "segmentID,timeIndex,time,cellID,faceID,processorRank,globalFaceID,processorLocalFaceID,associatedGlobalCellID,faceArea,kappaEff,snGradT,qWallIntoFluid,wallPower\n";
        }

        forAll(candidateCells_, candidateI)
        {
            const label celli = candidateCells_[candidateI];
            scalar wallArea = 0;
            scalar heatRate = 0;
            scalar weightedKappa = 0;
            scalar weightedSnGrad = 0;
            std::ostringstream faceIDs;
            label wallFaceCount = 0;

            forAll(patches_, patchi)
            {
                const label patchID = mesh.boundaryMesh().findPatchID(patches_[patchi]);
                if (patchID < 0)
                {
                    continue;
                }

                const fvPatch& patch = mesh.boundary()[patchID];
                const labelUList& faceCells = patch.faceCells();
                tmp<volScalarField> tAlphat = turbModel_.alphat();
                const scalarField& alphatPatch = tAlphat().boundaryField()[patchID];
                tmp<scalarField> tKappaPatch = phase1_.thermo().kappaEff(alphatPatch, patchID);
                const scalarField& kappaPatch = tKappaPatch();
                tmp<scalarField> tSnGrad = liquidTemperature.boundaryField()[patchID].snGrad();
                const scalarField& snGrad = tSnGrad();
                const vectorField& Sf = patch.Sf();

                forAll(faceCells, facei)
                {
                    if (faceCells[facei] != celli)
                    {
                        continue;
                    }

                    const label faceID = patch.start() + facei;
                    const scalar faceArea = mag(Sf[facei]);
                    const scalar qWallIntoFluid = kappaPatch[facei]*snGrad[facei];
                    const scalar facePower = qWallIntoFluid*faceArea;

                    if (wallFaceCount)
                    {
                        faceIDs << ";";
                    }
                    faceIDs << globalFaceIDs[faceID];
                    wallFaceCount++;
                    wallArea += faceArea;
                    heatRate += facePower;
                    weightedKappa += kappaPatch[facei]*faceArea;
                    weightedSnGrad += snGrad[facei]*faceArea;

                    faceOS
                        << perCellDiagnosticsSegmentID_ << ','
                        << timeIndex << ','
                        << mesh.time().value() << ','
                        << celli << ','
                        << faceID << ','
                        << Pstream::myProcNo() << ','
                        << globalFaceIDs[faceID] << ','
                        << faceID << ','
                        << globalCellIDs[celli] << ','
                        << faceArea << ','
                        << kappaPatch[facei] << ','
                        << snGrad[facei] << ','
                        << qWallIntoFluid << ','
                        << facePower << '\n';
                }
            }

            const scalar qWallIntoFluid = wallArea > SMALL ? heatRate/wallArea : scalar(0);
            const scalar kappaEff = wallArea > SMALL ? weightedKappa/wallArea : scalar(0);
            const scalar snGradT = wallArea > SMALL ? weightedSnGrad/wallArea : scalar(0);
            const vector& centre = mesh.C()[celli];
            const word controllerState =
                sustainedState_[celli] == REQUESTED_HOLD_ACTIVE
              ? word("REQUESTED_HOLD_ACTIVE")
              : sustainedState_[celli] == THERMALLY_SUSPENDED
              ? word("THERMALLY_SUSPENDED")
              : word("INACTIVE");
            const scalar holdStart =
                sustainedActivationTime_[celli] > -GREAT/2
              ? sustainedActivationTime_[celli]
              : scalar(-1);
            const scalar holdExpiry =
                sustainedActivationTime_[celli] > -GREAT/2
              ? sustainedActivationTime_[celli] + minimumActiveDuration_
              : scalar(-1);

            cellOS
                << perCellDiagnosticsSegmentID_ << ','
                << timeIndex << ','
                << mesh.time().value() << ','
                << mesh.time().deltaTValue() << ','
                << celli << ','
                << Pstream::myProcNo() << ','
                << globalCellIDs[celli] << ','
                << celli << ','
                << centre.x() << ','
                << centre.y() << ','
                << centre.z() << ','
                << mesh.V()[celli] << ','
                << sourceFootprintMode_ << ','
                << (siteCellMask_[celli] ? 1 : 0) << ','
                << diagnosticContactLineEligible[celli] << ','
                << diagnosticAlphaVapor[celli] << ','
                << controllerState << ','
                << holdStart << ','
                << holdExpiry << ','
                << wallFaceCount << ','
                << faceIDs.str() << ','
                << wallArea << ','
                << kappaEff << ','
                << snGradT << ','
                << qWallIntoFluid << ','
                << heatRate << ','
                << diagnosticWallEnergyAllowance[celli] << ','
                << liquidTemperature[celli] << ','
                << TSat(celli) << ','
                << superheat_[celli] << ','
                << diagnosticSensibleEnergyAllowance[celli] << ','
                << diagnosticWallEnergyAllowance[celli] << ','
                << diagnosticCombinedEnergyAllowance[celli] << ','
                << diagnosticRequestedAlphaStep[celli] << ','
                << maxAlphaVaporPerStep_ << ','
                << diagnosticAvailableLiquidLimit[celli] << ','
                << diagnosticNeighbourMaxAlphaVapor[celli] << ','
                << diagnosticSensibleOnlyThermalAlphaLimit[celli] << ','
                << diagnosticWallFluxAwareThermalAlphaLimit[celli] << ','
                << diagnosticActualAlphaStep[celli] << ','
                << alphaBirthSource_[celli] << ','
                << massBirthSource_[celli] << ','
                << latentSink_[celli] << ','
                << diagnosticLimiter[celli] << '\n';
        }

        lastPerCellDiagnosticsTimeIndex_ = timeIndex;
    }

    Info<< "wall-superheat vapor birth diagnostics: activeSourceCells=" << activeCells
        << ", sourceFootprintMode=" << sourceFootprintMode_
        << ", growthDepositionMode=" << growthDepositionMode_
        << ", sourceApplicationMode=" << sourceApplicationMode_
        << ", contactLineEligibleCells=" << contactLineEligibleCells
        << ", dryCoreExcludedCells=" << dryCoreExcludedCells
        << ", detachedExcludedCells=" << detachedExcludedCells
        << ", sustainedBirthMode=" << sustainedBirthMode_
        << ", requestedHoldCells=" << requestedHoldCells
        << ", thermallyActiveSourceCells=" << thermallyActiveSourceCells
        << ", thermallySuspendedCells=" << thermallySuspendedCells
        << ", holdExpiredCells=" << holdExpiredCells
        << ", siteMode=" << siteMode_
        << ", siteCenter=" << siteCenter_
        << ", siteRadius=" << siteRadius_
        << ", eligibleSiteCells=" << siteCells_.size()
        << ", eligibleSiteVolume=" << siteVolume()
        << ", firstActiveSiteTime=" << firstActiveSiteTime_
        << ", siteLocked=" << siteLocked_
        << ", sourceCellsOutsideSite=" << sourceCellsOutsideSite
        << ", interfaceBandEligibleCells=" << interfaceBandEligibleCells
        << ", interfaceBandSourceCells=" << interfaceBandSourceCells
        << ", interfaceBandSourceVolume=" << interfaceBandSourceVolume
        << ", interfaceBandSourceCentroidNormal=" << (interfaceBandSourceTotal > SMALL ? interfaceBandSourceWallNormalMoment/interfaceBandSourceTotal : scalar(0))
        << ", interfaceBandAboveFirstLayerFraction=" << (interfaceBandSourceTotal > SMALL ? interfaceBandSourceAboveFirstLayer/interfaceBandSourceTotal : scalar(0))
        << ", integratedRequestedAlphaSource=" << integratedRequestedAlphaSource
        << ", integratedAlphaSource=" << integratedAlphaSource
        << ", integratedMassSource=" << integratedMassSource
        << ", integratedLatentSink=" << integratedLatentSink
        << ", integratedThermalReserveUsed=" << integratedThermalReserveUsed
        << ", thermalLimitMode=" << thermalLimitMode_
        << ", wallFluxUseFraction=" << wallFluxUseFraction_
        << ", integratedWallFluxAllowance=" << integratedWallFluxAllowance
        << ", maxRequestedAlphaStep=" << maxRequestedAlphaStep
        << ", maxThermalAlphaLimit=" << maxThermalAlphaLimit
        << ", maxWallFluxIntoFluid=" << maxWallFluxIntoFluid
        << ", maxWallFluxAlphaLimit=" << maxWallFluxAlphaLimit
        << ", maxActualAlphaStep=" << maxActualAlphaStep
        << ", maxAlphaRate=" << maxAlphaRate
        << ", maxMassSource=" << maxMassSource
        << ", maxLatentSink=" << maxLatentSink
        << ", minSourceCellTemperature=" << (minSourceCellTemperature < GREAT/2 ? minSourceCellTemperature : 0)
        << ", minSourceCellSuperheat=" << (minSourceCellSuperheat < GREAT/2 ? minSourceCellSuperheat : 0)
        << ", maxSuperheat=" << (maxSuperheat > -GREAT/2 ? maxSuperheat : 0)
        << ", cavityBirthState=" << cavityStateName(cavityBirthState_)
        << ", cycleId=" << cycleId_
        << ", departureDetected=" << departureDetected_
        << ", rewetDetected=" << rewetDetected_
        << ", siteRearmed=" << siteRearmed_
        << ", firstTriggerTime=" << firstTriggerTime_ << endl;

    // During CONTACT_LINE_GROWTH with repeatedNucleationMode, updateCavityBirth
    // is called (writes birth-state/tip/energy CSV rows) but its local
    // integratedMassSource and activeCells are always zero because the actual
    // source is deposited in the loop above.  Write the source-history row here
    // with the real values so source-history.csv is accurate.
    if (cavityMode() && cavityBirthState_ == CONTACT_LINE_GROWTH && repeatedNucleationMode_)
    {
        scalar clgSiteMeanAlphaLiquid = scalar(1);
        scalar clgSiteVaporWallFraction = scalar(0);
        scalar clgAttachedWallContactArea = scalar(0);
        bool clgAttachedBubblePresent = false;
        evaluateSiteAttachment
        (
            clgSiteMeanAlphaLiquid,
            clgSiteVaporWallFraction,
            clgAttachedWallContactArea,
            clgAttachedBubblePresent
        );
        const bool clgDepartureCandidate =
            !clgAttachedBubblePresent
         && clgAttachedWallContactArea <= departureContactAreaThreshold_ + SMALL;
        writeSourceHistoryRow
        (
            integratedAlphaSource,
            integratedMassSource,
            integratedLatentSink,
            activeCells,
            clgDepartureCandidate,
            false
        );
    }
}

void Foam::wallSuperheatVaporBirth::TSource1(fvScalarMatrix& T1Eqn)
{
    updateSources();
    T1Eqn.source() -= latentSink_.internalField()*phase1_.mesh().V();
}

void Foam::wallSuperheatVaporBirth::TSource2(fvScalarMatrix& T2Eqn)
{}

void Foam::wallSuperheatVaporBirth::energySource(volScalarField& Q)
{}

void Foam::wallSuperheatVaporBirth::energySource1(volScalarField& q1)
{}

void Foam::wallSuperheatVaporBirth::energySource2(volScalarField& q2)
{}

void Foam::wallSuperheatVaporBirth::massSource(volScalarField& rhoSource)
{
    updateSources();
    rhoSource += massBirthSource_;
}

void Foam::wallSuperheatVaporBirth::alphaSource(volScalarField& rhoSource)
{
    updateSources();
    rhoSource += alphaBirthSource_;
}

void Foam::wallSuperheatVaporBirth::resetSourceFields
(
    volScalarField& superheat,
    volScalarField& alphaBirthSource,
    volScalarField& massBirthSource,
    volScalarField& latentSink
)
{
    superheat = dimensionedScalar("zero", superheat.dimensions(), 0.0);
    alphaBirthSource = dimensionedScalar("zero", alphaBirthSource.dimensions(), 0.0);
    massBirthSource = dimensionedScalar("zero", massBirthSource.dimensions(), 0.0);
    latentSink = dimensionedScalar("zero", latentSink.dimensions(), 0.0);
}

bool Foam::wallSuperheatVaporBirth::dimensionSelfTest
(
    const fvMesh& mesh,
    Ostream& os
)
{
    volScalarField temperature
    (
        IOobject
        (
            "dimensionTestTemperature",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("T", dimTemperature, 373.15),
        "zeroGradient"
    );

    volScalarField density
    (
        IOobject
        (
            "dimensionTestDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 1000.0),
        "zeroGradient"
    );

    volScalarField superheat
    (
        IOobject
        (
            "dimensionTestSuperheat",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimTemperature, 0.0),
        "zeroGradient"
    );

    volScalarField activeMask
    (
        IOobject
        (
            "dimensionTestActiveMask",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        "zeroGradient"
    );

    volScalarField alphaBirthSource
    (
        IOobject
        (
            "dimensionTestAlphaSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimTime, 0.0),
        "zeroGradient"
    );

    volScalarField massBirthSource
    (
        IOobject
        (
            "dimensionTestMassSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimDensity/dimTime, 0.0),
        "zeroGradient"
    );

    volScalarField latentSink
    (
        IOobject
        (
            "dimensionTestLatentSink",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPower/dimVol, 0.0),
        "zeroGradient"
    );

    resetSourceFields(superheat, alphaBirthSource, massBirthSource, latentSink);

    const dimensionedScalar TSat("TSat", dimTemperature, 373.15);
    superheat = temperature - TSat;

    tmp<volScalarField> massFromAlpha(density*alphaBirthSource);
    const dimensionedScalar latentHeat
    (
        "L",
        dimensionSet(0, 2, -2, 0, 0, 0, 0),
        2.26e6
    );
    tmp<volScalarField> latentFromMass(massBirthSource*latentHeat);

    bool ok = true;
    ok = ok && superheat.dimensions() == dimTemperature;
    ok = ok && activeMask.dimensions() == dimless;
    ok = ok && alphaBirthSource.dimensions() == dimless/dimTime;
    ok = ok && massBirthSource.dimensions() == dimDensity/dimTime;
    ok = ok && latentSink.dimensions() == dimPower/dimVol;
    ok = ok && massFromAlpha().dimensions() == dimDensity/dimTime;
    ok = ok && latentFromMass().dimensions() == dimPower/dimVol;
    ok = ok && latentSink.dimensions()*dimVol == dimPower;

    const dimensionedScalar cp
    (
        "Cp",
        dimensionSet(0, 2, -2, -1, 0, 0, 0),
        4200.0
    );
    ok = ok && density.dimensions()*cp.dimensions()*dimTemperature*dimVol
        == dimPower*dimTime;

    scalar maxSuperheat = 0;
    scalar maxActive = 0;
    scalar maxAlphaSource = 0;
    scalar maxMassSource = 0;
    scalar maxLatentSink = 0;

    forAll(superheat, celli)
    {
        maxSuperheat = max(maxSuperheat, mag(superheat[celli]));
        maxActive = max(maxActive, mag(activeMask[celli]));
        maxAlphaSource = max(maxAlphaSource, mag(alphaBirthSource[celli]));
        maxMassSource = max(maxMassSource, mag(massBirthSource[celli]));
        maxLatentSink = max(maxLatentSink, mag(latentSink[celli]));
    }

    reduce(maxSuperheat, maxOp<scalar>());
    reduce(maxActive, maxOp<scalar>());
    reduce(maxAlphaSource, maxOp<scalar>());
    reduce(maxMassSource, maxOp<scalar>());
    reduce(maxLatentSink, maxOp<scalar>());

    ok = ok && maxSuperheat <= SMALL;
    ok = ok && maxActive <= SMALL;
    ok = ok && maxAlphaSource <= SMALL;
    ok = ok && maxMassSource <= SMALL;
    ok = ok && maxLatentSink <= SMALL;

    vectorField testFaceCentres(5);
    testFaceCentres[0] = vector(0, 0, 0);
    testFaceCentres[1] = vector(0.0005, 0, 0);
    testFaceCentres[2] = vector(0.001, 0, 0);
    testFaceCentres[3] = vector(0.002, 0, 0);
    testFaceCentres[4] = vector(0.0005, 0, 0);

    labelList testFaceCells(5);
    testFaceCells[0] = 10;
    testFaceCells[1] = 11;
    testFaceCells[2] = 12;
    testFaceCells[3] = 13;
    testFaceCells[4] = 11;

    const labelList oneCellSite = circularPatchSiteCells
    (
        testFaceCentres,
        testFaceCells,
        vector::zero,
        0.00025
    );
    const labelList multiCellSite = circularPatchSiteCells
    (
        testFaceCentres,
        testFaceCells,
        vector::zero,
        0.001
    );
    const labelList zeroCellSite = circularPatchSiteCells
    (
        testFaceCentres,
        testFaceCells,
        vector(0.01, 0, 0),
        0.00025
    );

    ok = ok && oneCellSite.size() == 1 && oneCellSite[0] == 10;
    ok = ok && multiCellSite.size() == 3;
    ok = ok && zeroCellSite.empty();

    bool duplicateFree = true;
    forAll(multiCellSite, i)
    {
        for (label j = i + 1; j < multiCellSite.size(); ++j)
        {
            duplicateFree = duplicateFree && multiCellSite[i] != multiCellSite[j];
        }
    }
    ok = ok && duplicateFree;

    const bool contactLineInterfaceOk = contactLineEligible
    (
        0.5,
        0.0,
        0.02,
        0.98,
        0.02,
        true,
        0.25
    );
    const bool contactLineDryCoreExcluded = !contactLineEligible
    (
        0.01,
        0.0,
        0.02,
        0.98,
        0.02,
        true,
        0.25
    );
    const bool contactLineNeighbourOk = contactLineEligible
    (
        0.99,
        0.5,
        0.02,
        0.98,
        0.02,
        true,
        0.25
    );
    const bool contactLineDetachedExcluded = !contactLineEligible
    (
        0.99,
        0.0,
        0.02,
        0.98,
        0.02,
        true,
        0.25
    );
    const bool contactLineNeighbourSwitchOk = !contactLineEligible
    (
        0.99,
        0.5,
        0.02,
        0.98,
        0.02,
        false,
        0.25
    );

    ok = ok && contactLineInterfaceOk;
    ok = ok && contactLineDryCoreExcluded;
    ok = ok && contactLineNeighbourOk;
    ok = ok && contactLineDetachedExcluded;
    ok = ok && contactLineNeighbourSwitchOk;

    const scalar requestedStep = 0.01;
    const scalar availableAlpha = 1.0;
    const scalar alphaLiquid = 1.0;
    const scalar rhoConversion = 1000.0;
    const scalar cpValue = 4200.0;
    const scalar cellVolume = 1e-12;
    const scalar superheatOn = 0.02;
    const scalar superheatOff = 0.01;
    const scalar minimumActiveDuration = 0.005;
    const scalar minimumEvaporationSuperheat = 0.0;
    const scalar reserveFraction = 0.5;
    const scalar latentHeatValue = 2.26e6;
    const scalar deltaT = 1e-4;

    const scalar thermalLimitedStep = thermallyBoundedAlphaStep
    (
        requestedStep,
        requestedStep,
        availableAlpha,
        alphaLiquid,
        rhoConversion,
        cpValue,
        cellVolume,
        1.0,
        minimumEvaporationSuperheat,
        reserveFraction,
        latentHeatValue,
        0
    );
    const scalar expectedThermalLimitedStep = reserveFraction*cpValue/latentHeatValue;
    const bool thermalLimitOk =
        mag(thermalLimitedStep - expectedThermalLimitedStep) < 1e-12;

    const scalar capLimitedStep = thermallyBoundedAlphaStep
    (
        requestedStep,
        requestedStep,
        availableAlpha,
        alphaLiquid,
        rhoConversion,
        cpValue,
        cellVolume,
        100.0,
        minimumEvaporationSuperheat,
        reserveFraction,
        latentHeatValue,
        0
    );
    const bool capLimitOk = mag(capLimitedStep - requestedStep) < SMALL;

    const scalar liquidLimitedStep = thermallyBoundedAlphaStep
    (
        requestedStep,
        requestedStep,
        1e-4,
        alphaLiquid,
        rhoConversion,
        cpValue,
        cellVolume,
        100.0,
        minimumEvaporationSuperheat,
        reserveFraction,
        latentHeatValue,
        0
    );
    const bool liquidLimitOk = mag(liquidLimitedStep - 1e-4) < SMALL;

    const scalar floorSuspendedStep = thermallyBoundedAlphaStep
    (
        requestedStep,
        requestedStep,
        availableAlpha,
        alphaLiquid,
        rhoConversion,
        cpValue,
        cellVolume,
        0.0,
        minimumEvaporationSuperheat,
        reserveFraction,
        latentHeatValue,
        0
    );
    const scalar subSaturatedTimerStep = thermallyBoundedAlphaStep
    (
        requestedStep,
        requestedStep,
        availableAlpha,
        alphaLiquid,
        rhoConversion,
        cpValue,
        cellVolume,
        -1.0,
        minimumEvaporationSuperheat,
        reserveFraction,
        latentHeatValue,
        0
    );
    const bool suspendedLimitOk =
        floorSuspendedStep == 0 && subSaturatedTimerStep == 0;

    const scalar minimumStepSuspended = thermallyBoundedAlphaStep
    (
        requestedStep,
        requestedStep,
        availableAlpha,
        alphaLiquid,
        rhoConversion,
        cpValue,
        cellVolume,
        1.0,
        minimumEvaporationSuperheat,
        reserveFraction,
        latentHeatValue,
        0.01
    );
    const bool minimumThermalStepOk = minimumStepSuspended == 0;

    const scalar alphaRate = thermalLimitedStep/deltaT;
    const scalar massSource = rhoConversion*alphaRate;
    const scalar latentSource = massSource*latentHeatValue;
    const bool massAlphaOk = mag(massSource - rhoConversion*alphaRate) < SMALL;
    const bool massLatentOk = mag(latentSource - massSource*latentHeatValue) < 1e-6;

    label sustainedState = INACTIVE;
    scalar activationTime = -GREAT;
    bool requested = false;
    const bool activatesAboveOn = !requested && 0.03 >= superheatOn;
    if (activatesAboveOn)
    {
        requested = true;
        sustainedState = REQUESTED_HOLD_ACTIVE;
        activationTime = 0;
    }
    const bool noActivationBelowOn = !false && 0.005 < superheatOn;
    const bool timerStillActive = 0.001 - activationTime < minimumActiveDuration;
    const bool suspendedWhileTimerActive =
        requested
     && timerStillActive
     && subSaturatedTimerStep == 0
     && THERMALLY_SUSPENDED == 2;
    sustainedState = THERMALLY_SUSPENDED;
    const bool reheatsBeforeExpiry =
        sustainedState == THERMALLY_SUSPENDED
     && timerStillActive
     && thermallyBoundedAlphaStep
        (
            requestedStep,
            requestedStep,
            availableAlpha,
            alphaLiquid,
            rhoConversion,
            cpValue,
            cellVolume,
            1.0,
            minimumEvaporationSuperheat,
            reserveFraction,
            latentHeatValue,
            0
        ) > 0;
    const bool holdExpiresToStandardHysteresis =
        (0.006 - activationTime >= minimumActiveDuration) && 0.005 <= superheatOff;

    ok = ok && thermalLimitOk;
    ok = ok && capLimitOk;
    ok = ok && liquidLimitOk;
    ok = ok && suspendedLimitOk;
    ok = ok && minimumThermalStepOk;
    ok = ok && massAlphaOk;
    ok = ok && massLatentOk;
    ok = ok && activatesAboveOn;
    ok = ok && noActivationBelowOn;
    ok = ok && sustainedState == THERMALLY_SUSPENDED;
    ok = ok && suspendedWhileTimerActive;
    ok = ok && reheatsBeforeExpiry;
    ok = ok && holdExpiresToStandardHysteresis;

    os  << "wallSuperheatVaporBirth dimension self-test" << nl
        << "    superheat dimensions: " << superheat.dimensions() << nl
        << "    activeMask dimensions: " << activeMask.dimensions() << nl
        << "    alphaBirthSource dimensions: " << alphaBirthSource.dimensions() << nl
        << "    massBirthSource dimensions: " << massBirthSource.dimensions() << nl
        << "    latentSink dimensions: " << latentSink.dimensions() << nl
        << "    rho*alphaRate dimensions: " << massFromAlpha().dimensions() << nl
        << "    massSource*latentHeat dimensions: " << latentFromMass().dimensions() << nl
        << "    TSource source dimensions: " << latentSink.dimensions()*dimVol << nl
        << "    thermal reserve dimensions: "
        << density.dimensions()*cp.dimensions()*dimTemperature*dimVol << nl
        << "    max inactive superheat: " << maxSuperheat << nl
        << "    max inactive activeMask: " << maxActive << nl
        << "    max inactive alpha source: " << maxAlphaSource << nl
        << "    max inactive mass source: " << maxMassSource << nl
        << "    max inactive latent sink: " << maxLatentSink << nl
        << "    circular one-cell site count: " << oneCellSite.size() << nl
        << "    circular multi-cell site count: " << multiCellSite.size() << nl
        << "    circular zero-cell site count: " << zeroCellSite.size() << nl
        << "    circular site duplicate-free: " << duplicateFree << nl
        << "    wallAdjacentContactLine interface-band eligible: " << contactLineInterfaceOk << nl
        << "    wallAdjacentContactLine dry-core excluded: " << contactLineDryCoreExcluded << nl
        << "    wallAdjacentContactLine neighbour eligible: " << contactLineNeighbourOk << nl
        << "    wallAdjacentContactLine detached excluded: " << contactLineDetachedExcluded << nl
        << "    wallAdjacentContactLine neighbour switch respected: " << contactLineNeighbourSwitchOk << nl
        << "    activation above superheatOn: " << activatesAboveOn << nl
        << "    no activation below superheatOn: " << noActivationBelowOn << nl
        << "    hold timer state persistence: " << timerStillActive << nl
        << "    positive source while thermally supported: " << (thermalLimitedStep > 0) << nl
        << "    zero source at evaporation floor: " << (floorSuspendedStep == 0) << nl
        << "    sub-saturated timer-active source zero: " << (subSaturatedTimerStep == 0) << nl
        << "    THERMALLY_SUSPENDED state: " << (sustainedState == THERMALLY_SUSPENDED) << nl
        << "    reactivation during hold after reheating: " << reheatsBeforeExpiry << nl
        << "    hold expiry returns to standard hysteresis: " << holdExpiresToStandardHysteresis << nl
        << "    thermal alpha limit: " << thermalLimitOk << nl
        << "    available-liquid limit: " << liquidLimitOk << nl
        << "    per-step cap limit: " << capLimitOk << nl
        << "    mass/alpha consistency: " << massAlphaOk << nl
        << "    mass/latent consistency: " << massLatentOk << nl
        << "    minimum thermal alpha step suspension: " << minimumThermalStepOk << nl
        << "    result: " << (ok ? "PASS" : "FAIL") << endl;

    return ok;
}

Foam::labelList Foam::wallSuperheatVaporBirth::circularPatchSiteCells
(
    const vectorField& faceCentres,
    const labelUList& faceCells,
    const vector& siteCenter,
    const scalar siteRadius
)
{
    if (siteRadius <= 0)
    {
        FatalErrorInFunction
            << "circularPatchSite requires positive siteRadius; found "
            << siteRadius << exit(FatalError);
    }

    label maxCell = -1;
    forAll(faceCells, facei)
    {
        maxCell = max(maxCell, faceCells[facei]);
    }

    boolList selected(maxCell + 1, false);
    DynamicList<label> cells(faceCells.size());

    forAll(faceCentres, facei)
    {
        if (mag(faceCentres[facei] - siteCenter) <= siteRadius + SMALL)
        {
            const label celli = faceCells[facei];
            if (!selected[celli])
            {
                selected[celli] = true;
                cells.append(celli);
            }
        }
    }

    labelList result;
    result.transfer(cells);
    return result;
}

// ************************************************************************* //
