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
#include "mathematicalConstants.H"

#include <cmath>
#include <cstdint>

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
    maxPreexistingHandoffVaporFraction_
    (
        modelDict().lookupOrDefault<scalar>
        (
            "maxPreexistingHandoffVaporFraction",
            1.0
        )
    ),
    maximumHandoffDuration_
    (
        modelDict().lookupOrDefault<scalar>("maximumHandoffDuration", -1.0)
    ),
    creationSteps_(modelDict().lookupOrDefault<label>("creationSteps", 30)),
    maximumCreationSteps_
    (
        modelDict().lookupOrDefault<label>("maximumCreationSteps", 300)
    ),
    creationDuration_
    (
        modelDict().lookupOrDefault<scalar>("creationDuration", -1.0)
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
    wallBoilingClosure_
    (
        modelDict().lookupOrDefault<word>("wallBoilingClosure", "none")
    ),
    wallBoilingCoefficient_
    (
        modelDict().lookupOrDefault<scalar>("wallBoilingCoefficient", 0.0)
    ),
    wallBoilingExponent_
    (
        modelDict().lookupOrDefault<scalar>("wallBoilingExponent", 3.0)
    ),
    wallBoilingMaxFlux_
    (
        modelDict().lookupOrDefault<scalar>("wallBoilingMaxFlux", -1.0)
    ),
    wallBoilingFluxMode_
    (
        modelDict().lookupOrDefault<word>("wallBoilingFluxMode", "total")
    ),
    wallBoilingFluxFieldName_
    (
        modelDict().lookupOrDefault<word>
        (
            "wallBoilingFluxField",
            "qWallBoiling"
        )
    ),
    displacementShellCells_
    (
        modelDict().lookupOrDefault<scalar>("displacementShellCells", 4.0)
    ),
    departureDiameterDistribution_
    (
        modelDict().lookupOrDefault<word>
        (
            "departureDiameterDistribution",
            "fixed"
        )
    ),
    departureDiameterMedian_
    (
        modelDict().lookupOrDefault<scalar>
        (
            "departureDiameterMedian",
            2.0*departureRadius_
        )
    ),
    departureDiameterSigmaLn_
    (
        modelDict().lookupOrDefault<scalar>("departureDiameterSigmaLn", 0.0)
    ),
    departureDiameterMin_
    (
        modelDict().lookupOrDefault<scalar>
        (
            "departureDiameterMin",
            2.0*departureRadius_
        )
    ),
    departureDiameterMax_
    (
        modelDict().lookupOrDefault<scalar>
        (
            "departureDiameterMax",
            2.0*departureRadius_
        )
    ),
    randomSeed_(modelDict().lookupOrDefault<label>("randomSeed", 0)),
    rearmDelayScaling_
    (
        modelDict().lookupOrDefault<word>("rearmDelayScaling", "fixed")
    ),
    waitingTimeMean_
    (
        modelDict().lookupOrDefault<scalar>("waitingTimeMean", 0.0)
    ),
    wallHeatPartition_
    (
        modelDict().lookupOrDefault<word>("wallHeatPartition", "vapor")
    ),
    evaporationFraction_
    (
        modelDict().lookupOrDefault<scalar>("evaporationFraction", 1.0)
    ),
    evaporationFractionSuperheat_(),
    evaporationFractionValue_(),
    inactiveAreaHeat_
    (
        modelDict().lookupOrDefault<word>("inactiveAreaHeat", "none")
    ),
    writeDiagnostics_
    (
        modelDict().lookupOrDefault<Switch>("writeDiagnostics", true)
    ),
    writeHandoffAudit_
    (
        modelDict().lookupOrDefault<Switch>("writeHandoffAudit", false)
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
    displacementCells_(),
    captureArea_(),
    siteDepartureRadius_(),
    targetVaporVolume_(),
    targetVaporMass_(),
    targetLatentEnergy_(),
    siteState_(),
    activationDuration_(),
    storedEnergy_(),
    releaseTime_(),
    creationStartTime_(),
    initialSensibleCredit_(),
    chargedWallEnergy_(),
    nextActivationTime_(),
    creationStep_(),
    cycleId_(),
    createdVaporVolume_(),
    createdVaporMass_(),
    consumedLatentEnergy_(),
    boilingEnergy_(),
    cumulativeBoilingEnergy_(),
    cumulativeBoilingBubbleEnergy_(),
    cumulativeLiquidHeatingEnergy_(),
    cumulativeDisplacedLiquidVolume_(0.0),
    cumulativeBirthVaporVolume_(0.0),
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
    ),
    wallBoilingFlux_
    (
        IOobject
        (
            wallBoilingFluxFieldName_,
            phase1.mesh().time().timeName(),
            phase1.mesh(),
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        phase1.mesh(),
        dimensionedScalar("zero", dimPower/dimArea, 0.0),
        "calculated"
    ),
    wallLiquidHeat_
    (
        IOobject
        (
            word(diagnosticPrefix_ + "WallLiquidHeat"),
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
    if (wallHeatPartition_ == "table")
    {
        evaporationFractionSuperheat_ =
            scalarField(modelDict().lookup("evaporationFractionSuperheat"));
        evaporationFractionValue_ =
            scalarField(modelDict().lookup("evaporationFractionValue"));
    }

    validateControls();
    buildSiteStencils();

    const label nSites = siteCentres_.size();
    siteState_.setSize(nSites, WAITING_FOR_ACTIVATION);
    activationDuration_.setSize(nSites, 0.0);
    storedEnergy_.setSize(nSites, 0.0);
    releaseTime_.setSize(nSites, -GREAT);
    creationStartTime_.setSize(nSites, -GREAT);
    initialSensibleCredit_.setSize(nSites, 0.0);
    chargedWallEnergy_.setSize(nSites, 0.0);
    nextActivationTime_.setSize(nSites, 0.0);
    creationStep_.setSize(nSites, 0);
    cycleId_.setSize(nSites, 0);
    createdVaporVolume_.setSize(nSites, 0.0);
    createdVaporMass_.setSize(nSites, 0.0);
    consumedLatentEnergy_.setSize(nSites, 0.0);
    boilingEnergy_.setSize(nSites, 0.0);
    cumulativeBoilingEnergy_.setSize(nSites, 0.0);
    cumulativeBoilingBubbleEnergy_.setSize(nSites, 0.0);
    cumulativeLiquidHeatingEnergy_.setSize(nSites, 0.0);

    forAll(nextActivationTime_, siteI)
    {
        // Cycle 0 is the wait before a site's first activation; cycle 1 is
        // its first release (cycleId_ is incremented on activation).
        nextActivationTime_[siteI] =
            activationStagger_*siteI + drawWaitingTime(siteI, 0);
        if (departureDiameterDistribution_ != "fixed")
        {
            buildReleaseStencil(siteI, 0.5*drawDepartureDiameter(siteI, 1));
        }
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
            << ", maxPreexistingHandoffVaporFraction="
            << maxPreexistingHandoffVaporFraction_
            << ", maximumHandoffDuration=" << maximumHandoffDuration_
            << ", writeHandoffAudit=" << writeHandoffAudit_
            << ", continuitySourceMode=" << continuitySourceMode_
            << ", parallel=" << Pstream::parRun()
            << endl;
        Info<< "populationDetachedBubbleBirth wall boiling closure: "
            << wallBoilingClosure_
            << ", coefficient=" << wallBoilingCoefficient_
            << ", exponent=" << wallBoilingExponent_
            << ", maxFlux=" << wallBoilingMaxFlux_
            << ", fluxMode=" << wallBoilingFluxMode_
            << ", fluxField=" << wallBoilingFluxFieldName_
            << ", displacementShellCells=" << displacementShellCells_
            << ", creationDuration=" << creationDuration_
            << endl;
        Info<< "populationDetachedBubbleBirth heater-dependent controls: "
            << "departureDiameterDistribution=" << departureDiameterDistribution_
            << ", median=" << departureDiameterMedian_
            << ", sigmaLn=" << departureDiameterSigmaLn_
            << ", min=" << departureDiameterMin_
            << ", max=" << departureDiameterMax_
            << ", randomSeed=" << randomSeed_
            << ", rearmDelayScaling=" << rearmDelayScaling_
            << ", waitingTimeMean=" << waitingTimeMean_
            << ", wallHeatPartition=" << wallHeatPartition_
            << ", evaporationFraction=" << evaporationFraction_
            << ", evaporationFractionTablePoints="
            << evaporationFractionSuperheat_.size()
            << ", inactiveAreaHeat=" << inactiveAreaHeat_
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
    if
    (
        maxPreexistingHandoffVaporFraction_ < 0
     || maxPreexistingHandoffVaporFraction_ > 1
    )
    {
        FatalErrorInFunction
            << "maxPreexistingHandoffVaporFraction must be in [0,1]"
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
     && continuitySourceMode_ != "displacedLiquid"
     && continuitySourceMode_ != "none"
    )
    {
        FatalErrorInFunction
            << "continuitySourceMode must be vaporMass, displacedLiquid or none"
            << exit(FatalError);
    }
    if (creationDuration_ != -1.0 && creationDuration_ <= 0)
    {
        FatalErrorInFunction
            << "creationDuration must be positive or -1 (use creationSteps)"
            << exit(FatalError);
    }
    if (displacementShellCells_ <= 0)
    {
        FatalErrorInFunction
            << "displacementShellCells must be positive"
            << exit(FatalError);
    }
    if
    (
        wallBoilingClosure_ != "none"
     && wallBoilingClosure_ != "powerLaw"
    )
    {
        FatalErrorInFunction
            << "wallBoilingClosure must be none or powerLaw; found "
            << wallBoilingClosure_
            << exit(FatalError);
    }
    if
    (
        wallBoilingClosure_ == "powerLaw"
     && (wallBoilingCoefficient_ <= 0 || wallBoilingExponent_ < 1)
    )
    {
        FatalErrorInFunction
            << "wallBoilingClosure powerLaw requires "
            << "wallBoilingCoefficient > 0 and wallBoilingExponent >= 1"
            << exit(FatalError);
    }
    if
    (
        wallBoilingFluxMode_ != "total"
     && wallBoilingFluxMode_ != "additional"
    )
    {
        FatalErrorInFunction
            << "wallBoilingFluxMode must be total or additional; found "
            << wallBoilingFluxMode_
            << exit(FatalError);
    }
    if
    (
        departureDiameterDistribution_ != "fixed"
     && departureDiameterDistribution_ != "logNormal"
    )
    {
        FatalErrorInFunction
            << "departureDiameterDistribution must be fixed or logNormal; "
            << "found " << departureDiameterDistribution_
            << exit(FatalError);
    }
    if
    (
        departureDiameterDistribution_ == "logNormal"
     && (
            departureDiameterMin_ <= 0
         || departureDiameterMax_ < departureDiameterMin_
         || departureDiameterMedian_ <= 0
         || departureDiameterSigmaLn_ < 0
        )
    )
    {
        FatalErrorInFunction
            << "logNormal departure diameters need 0 < departureDiameterMin"
            << " <= departureDiameterMax, departureDiameterMedian > 0 and"
            << " departureDiameterSigmaLn >= 0; found min="
            << departureDiameterMin_ << " max=" << departureDiameterMax_
            << " median=" << departureDiameterMedian_
            << " sigmaLn=" << departureDiameterSigmaLn_
            << exit(FatalError);
    }
    if (rearmDelayScaling_ != "fixed" && rearmDelayScaling_ != "diameter")
    {
        FatalErrorInFunction
            << "rearmDelayScaling must be fixed or diameter; found "
            << rearmDelayScaling_
            << exit(FatalError);
    }
    if (waitingTimeMean_ < 0)
    {
        FatalErrorInFunction
            << "waitingTimeMean must be non-negative; found "
            << waitingTimeMean_
            << exit(FatalError);
    }
    if
    (
        wallHeatPartition_ != "vapor"
     && wallHeatPartition_ != "constant"
     && wallHeatPartition_ != "table"
    )
    {
        FatalErrorInFunction
            << "wallHeatPartition must be vapor, constant or table; found "
            << wallHeatPartition_
            << exit(FatalError);
    }
    if (wallHeatPartition_ != "vapor" && wallBoilingClosure_ == "none")
    {
        FatalErrorInFunction
            << "wallHeatPartition " << wallHeatPartition_
            << " needs the wall boiling closure"
            << exit(FatalError);
    }
    if
    (
        inactiveAreaHeat_ != "none"
     && inactiveAreaHeat_ != "liquid"
     && inactiveAreaHeat_ != "liquidThenNearestSite"
     && inactiveAreaHeat_ != "liquidThenNearestSiteWithRoom"
    )
    {
        FatalErrorInFunction
            << "inactiveAreaHeat must be none, liquid, "
            << "liquidThenNearestSite or liquidThenNearestSiteWithRoom; found "
            << inactiveAreaHeat_
            << exit(FatalError);
    }
    if (inactiveAreaHeat_ != "none" && wallBoilingClosure_ == "none")
    {
        FatalErrorInFunction
            << "inactiveAreaHeat " << inactiveAreaHeat_
            << " needs the wall boiling closure"
            << exit(FatalError);
    }
    if
    (
        wallHeatPartition_ == "constant"
     && (evaporationFraction_ < 0 || evaporationFraction_ > 1)
    )
    {
        FatalErrorInFunction
            << "evaporationFraction must lie in [0, 1]; found "
            << evaporationFraction_
            << exit(FatalError);
    }
    if (wallHeatPartition_ == "table")
    {
        const label n = evaporationFractionSuperheat_.size();
        bool valid = n >= 2 && evaporationFractionValue_.size() == n;
        for (label i = 0; valid && i < n; ++i)
        {
            valid =
                evaporationFractionValue_[i] >= 0
             && evaporationFractionValue_[i] <= 1
             && (
                    i == 0
                 || evaporationFractionSuperheat_[i]
                  > evaporationFractionSuperheat_[i - 1]
                );
        }
        if (!valid)
        {
            FatalErrorInFunction
                << "wallHeatPartition table needs at least two points, "
                << "evaporationFractionSuperheat increasing and "
                << "evaporationFractionValue of the same length in [0, 1]"
                << exit(FatalError);
        }
    }
    if (maximumHandoffDuration_ != -1.0 && maximumHandoffDuration_ <= 0)
    {
        FatalErrorInFunction
            << "maximumHandoffDuration must be positive or -1 (disabled)"
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

    const label nSites = siteCentres_.size();
    siteDepartureRadius_.setSize(nSites, departureRadius_);
    siteFaceLabels_.setSize(nSites);
    sphereCells_.setSize(nSites);
    sphereTargetVapor_.setSize(nSites);
    shieldingCells_.setSize(nSites);
    displacementCells_.setSize(nSites);
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

        buildReleaseStencil(siteI, departureRadius_);

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
        label globalShieldCount = shieldingCells_[siteI].size();
        reduce(globalFaceCount, sumOp<label>());
        reduce(globalShieldCount, sumOp<label>());

        if
        (
            globalFaceCount == 0
         || globalShieldCount == 0
         || localCaptureArea <= SMALL
        )
        {
            FatalErrorInFunction
                << "Site " << siteI << " at " << site
                << " has an invalid capture or shielding stencil"
                << exit(FatalError);
        }
    }

    if (inactiveAreaHeat_ != "none")
    {
        // The option covers the faces sites own. Report how much of the
        // heater that is, so a gap (capture radii, or ownership centres
        // without a CFD site) shows up in the log.
        const scalar heaterArea = gSum(mag(faceAreaVectors));
        Info<< "POPULATION_INACTIVE_AREA_HEAT"
            << " mode=" << inactiveAreaHeat_
            << " siteArea=" << sum(captureArea_)
            << " heaterArea=" << heaterArea
            << endl;
    }
}

void Foam::populationDetachedBubbleBirth::buildReleaseStencil
(
    const label siteI,
    const scalar radius
)
{
    // The cells a site's bubble is created in, their target vapor fractions,
    // the latent heat that bubble needs, and the liquid shell that receives
    // the liquid it displaces, for a bubble of the given radius. Every
    // processor must call this for the same site at the same time: it ends
    // with global reductions.
    const fvMesh& mesh = phase1_.mesh();
    tmp<volScalarField> tRhoV = phase2_.thermo().rho();
    const volScalarField& rhoV = tRhoV();
    const vector& site = siteCentres_[siteI];
    siteDepartureRadius_[siteI] = radius;

    {
        const vector bubbleCentre =
            site + vector(0, 0, radius + detachmentGap_);
        DynamicList<label> localSphereCells;
        DynamicList<label> localDisplacementCells;
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
                radius - mag(mesh.C()[celli] - bubbleCentre);
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
            else if
            (
                continuitySourceMode_ == "displacedLiquid"
             && mag(mesh.C()[celli] - bubbleCentre)
             <= radius + displacementShellCells_*cellSize
            )
            {
                // Liquid just outside the new bubble. It receives the
                // volume of liquid the bubble pushes aside.
                localDisplacementCells.append(celli);
            }
        }
        sphereCells_[siteI].transfer(localSphereCells);
        displacementCells_[siteI].transfer(localDisplacementCells);
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

        label globalSphereCount = sphereCells_[siteI].size();
        reduce(globalSphereCount, sumOp<label>());
        label globalDisplacementCount = displacementCells_[siteI].size();
        reduce(globalDisplacementCount, sumOp<label>());
        if
        (
            continuitySourceMode_ == "displacedLiquid"
         && globalDisplacementCount == 0
        )
        {
            FatalErrorInFunction
                << "Site " << siteI << " at " << site
                << " has no cells around its bubble to receive displaced "
                << "liquid"
                << exit(FatalError);
        }

        if
        (
            globalSphereCount == 0
         || localVolume <= SMALL
         || localLatent <= SMALL
         || localCoreCells == 0
        )
        {
            FatalErrorInFunction
                << "Site " << siteI << " at " << site
                << " has an invalid release stencil for a bubble of radius "
                << radius << ": cells=" << globalSphereCount
                << " volume=" << localVolume
                << " coreCells=" << localCoreCells
                << exit(FatalError);
        }
    }
}

Foam::scalar Foam::populationDetachedBubbleBirth::hashUniform
(
    const label siteI,
    const label cycle,
    const label stream
) const
{
    // SplitMix64 over (seed, site, cycle, stream): a uniform number in
    // (0, 1) that every processor computes identically.
    auto mix = [](uint64_t x)
    {
        x += 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30))*0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27))*0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    };
    uint64_t h = mix(static_cast<uint64_t>(randomSeed_));
    h = mix(h ^ static_cast<uint64_t>(siteI));
    h = mix(h ^ static_cast<uint64_t>(cycle));
    h = mix(h ^ static_cast<uint64_t>(stream));
    // The top 53 bits, shifted off zero and one.
    return (static_cast<scalar>(h >> 11) + 1.0)/9007199254740994.0;
}

Foam::scalar Foam::populationDetachedBubbleBirth::drawDepartureDiameter
(
    const label siteI,
    const label cycle
) const
{
    if (departureDiameterDistribution_ == "fixed")
    {
        return 2.0*departureRadius_;
    }
    // Log-normal by Box-Muller, redrawn (on fresh streams) until the value
    // lies in [min, max].
    for (label attempt = 0; attempt < 64; ++attempt)
    {
        const scalar u1 = hashUniform(siteI, cycle, 2*attempt + 1);
        const scalar u2 = hashUniform(siteI, cycle, 2*attempt + 2);
        const scalar z =
            Foam::sqrt(-2.0*Foam::log(u1))
           *Foam::cos(constant::mathematical::twoPi*u2);
        const scalar d =
            departureDiameterMedian_*Foam::exp(departureDiameterSigmaLn_*z);
        if (d >= departureDiameterMin_ && d <= departureDiameterMax_)
        {
            return d;
        }
    }
    return min(max(departureDiameterMedian_, departureDiameterMin_), departureDiameterMax_);
}

Foam::scalar Foam::populationDetachedBubbleBirth::drawWaitingTime
(
    const label siteI,
    const label cycle
) const
{
    if (waitingTimeMean_ <= 0)
    {
        return 0;
    }
    return -waitingTimeMean_*Foam::log(hashUniform(siteI, cycle, 1000));
}

Foam::scalar Foam::populationDetachedBubbleBirth::rearmDelay
(
    const label siteI
) const
{
    if (rearmDelayScaling_ == "diameter")
    {
        return minimumRearmDelay_*siteDepartureRadius_[siteI]/departureRadius_;
    }
    return minimumRearmDelay_;
}

Foam::scalar Foam::populationDetachedBubbleBirth::siteEvaporationFraction
(
    const scalar superheat
) const
{
    if (wallHeatPartition_ == "constant")
    {
        return evaporationFraction_;
    }
    if (wallHeatPartition_ == "table")
    {
        const scalarField& x = evaporationFractionSuperheat_;
        const scalarField& y = evaporationFractionValue_;
        const label n = x.size();
        if (superheat <= x[0])
        {
            return y[0];
        }
        for (label i = 1; i < n; ++i)
        {
            if (superheat <= x[i])
            {
                const scalar w = (superheat - x[i - 1])/(x[i] - x[i - 1]);
                return (1 - w)*y[i - 1] + w*y[i];
            }
        }
        return y[n - 1];
    }
    return 1;
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
    creationStartTime_[siteI] = -GREAT;
    initialSensibleCredit_[siteI] = 0.0;
    chargedWallEnergy_[siteI] = 0.0;
    // The delay follows the bubble just released (or abandoned). The
    // release stencil is then rebuilt for the next cycle's bubble, so that
    // everything the site measures while it waits and charges refers to it.
    nextActivationTime_[siteI] =
        timeValue + rearmDelay(siteI) + drawWaitingTime(siteI, cycleId_[siteI]);
    if (departureDiameterDistribution_ != "fixed")
    {
        buildReleaseStencil
        (
            siteI,
            0.5*drawDepartureDiameter(siteI, cycleId_[siteI] + 1)
        );
    }
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
    wallLiquidHeat_ *= scalar(0);

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

    const bool boilingClosure = wallBoilingClosure_ != "none";
    fvPatchScalarField& boilingFluxPatch =
        wallBoilingFlux_.boundaryFieldRef()[patchID_];
    boilingFluxPatch = scalar(0);
    scalar totalBoilingPower = 0;
    scalar totalConductionPower = 0;
    scalar totalLiquidHeatingPower = 0;
    const bool partitioned = wallHeatPartition_ != "vapor";
    const bool diagnosticsStep =
        timeIndex % diagnosticsInterval_ == 0 || mesh.time().writeTime();

    // inactiveAreaHeat liquidThenNearestSite(WithRoom): the faces of inactive
    // sites are handled here, before the site loop, so that the heat their
    // wall-cell liquid cannot take reaches an active site in this same step.
    // The law is the same as in the site loop below.
    const label nSites = siteCentres_.size();
    const bool withRoom = inactiveAreaHeat_ == "liquidThenNearestSiteWithRoom";
    const bool nearestSiteOverflow =
        boilingClosure
     && (inactiveAreaHeat_ == "liquidThenNearestSite" || withRoom);
    const label nPrepass = nearestSiteOverflow ? nSites : 0;
    boolList prepassActive(nPrepass, false);
    scalarField inactiveLiquidPower(nPrepass, 0.0);
    scalarField inactiveConductionPower(nPrepass, 0.0);
    scalarField receivedPower(nPrepass, 0.0);
    // Overflow that found no active site with room and went to the nearest
    // one regardless (liquidThenNearestSiteWithRoom).
    scalar overflowWithoutRoomPower = 0;
    if (nearestSiteOverflow)
    {
        forAll(siteCentres_, siteI)
        {
            scalar superheatArea = 0;
            scalar siteArea = 0;
            const labelList& faces = siteFaceLabels_[siteI];
            forAll(faces, faceListI)
            {
                const label faceI = faces[faceListI];
                const scalar area = mag(faceAreaVectors[faceI]);
                superheatArea +=
                    max(wallTemperature[faceI] - TSat(faceCells[faceI]), scalar(0))
                   *area;
                siteArea += area;
            }
            reduce(superheatArea, sumOp<scalar>());
            reduce(siteArea, sumOp<scalar>());
            const scalar meanSuperheat =
                siteArea > SMALL ? superheatArea/siteArea : scalar(0);
            prepassActive[siteI] =
                meanSuperheat + SMALL >= siteActivationSuperheat_[siteI];
        }

        // liquidThenNearestSiteWithRoom: the energy each active site still
        // lacks for the bubble it is charging, creating or will charge next.
        // A site that already holds that much has no room, and heat sent to
        // it would only wait there. The room shrinks as this pass sends heat.
        // Every processor holds the same site state, so all of them choose
        // the same targets.
        scalarField room(withRoom ? nSites : 0, 0.0);
        forAll(room, siteI)
        {
            const scalar need =
                targetLatentEnergy_[siteI]
              - (
                    siteState_[siteI] == CREATING_DETACHED_BUBBLE
                  ? consumedLatentEnergy_[siteI]
                  : scalar(0)
                );
            room[siteI] = need - boilingEnergy_[siteI] - storedEnergy_[siteI];
        }

        forAll(siteCentres_, siteI)
        {
            if (prepassActive[siteI])
            {
                continue;
            }
            // Nearest active site in the heater plane; -1 if none is active.
            // With room: the nearest active site that has room, if any.
            label targetI = -1;
            label roomTargetI = -1;
            scalar nearestDistanceSqr = GREAT;
            scalar nearestRoomDistanceSqr = GREAT;
            forAll(siteCentres_, otherI)
            {
                if (!prepassActive[otherI])
                {
                    continue;
                }
                const vector d = siteCentres_[otherI] - siteCentres_[siteI];
                const scalar distanceSqr = sqr(d.x()) + sqr(d.y());
                if (distanceSqr < nearestDistanceSqr)
                {
                    nearestDistanceSqr = distanceSqr;
                    targetI = otherI;
                }
                if
                (
                    withRoom
                 && room[otherI] > 0
                 && distanceSqr < nearestRoomDistanceSqr
                )
                {
                    nearestRoomDistanceSqr = distanceSqr;
                    roomTargetI = otherI;
                }
            }

            scalar liquidPower = 0;
            scalar overflowPower = 0;
            scalar conductionPower = 0;
            const labelList& faces = siteFaceLabels_[siteI];
            forAll(faces, faceListI)
            {
                const label faceI = faces[faceListI];
                const label celli = faceCells[faceI];
                const scalar area = mag(faceAreaVectors[faceI]);
                const scalar superheat =
                    max(wallTemperature[faceI] - TSat(celli), scalar(0));
                const scalar conductionFlux =
                    max(kappaPatch[faceI]*snGrad[faceI], scalar(0));
                scalar lawFlux =
                    wallBoilingCoefficient_
                   *Foam::pow(superheat, wallBoilingExponent_);
                if (wallBoilingMaxFlux_ > 0)
                {
                    lawFlux = min(lawFlux, wallBoilingMaxFlux_);
                }
                scalar boilingFlux =
                    wallBoilingFluxMode_ == "total"
                  ? max(lawFlux - conductionFlux, scalar(0))
                  : lawFlux;
                boilingFlux *=
                    max(min(phase1_[celli], scalar(1)), scalar(0));
                conductionPower += conductionFlux*area;

                if (liquidTemperature[celli] < wallTemperature[faceI])
                {
                    // The wall-cell liquid takes the heat.
                    boilingFluxPatch[faceI] = -boilingFlux;
                    wallLiquidHeat_[celli] += boilingFlux*area/mesh.V()[celli];
                    liquidPower += boilingFlux*area;
                }
                else if (targetI >= 0)
                {
                    // The liquid is full: the heat funds the bubbles of an
                    // active site instead (the nearest, or with room the
                    // nearest that has room; chosen below).
                    boilingFluxPatch[faceI] = -boilingFlux;
                    overflowPower += boilingFlux*area;
                }
            }
            reduce(liquidPower, sumOp<scalar>());
            reduce(overflowPower, sumOp<scalar>());
            reduce(conductionPower, sumOp<scalar>());
            inactiveLiquidPower[siteI] = liquidPower;
            inactiveConductionPower[siteI] = conductionPower;
            if (withRoom && overflowPower > 0)
            {
                if (roomTargetI >= 0)
                {
                    targetI = roomTargetI;
                    room[targetI] -= overflowPower*deltaT;
                }
                else
                {
                    overflowWithoutRoomPower += overflowPower;
                }
            }
            if (targetI >= 0)
            {
                receivedPower[targetI] += overflowPower;
            }
        }
    }

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

        // Wall boiling closure. The flux leaves the solid through the
        // coupled boundary (qr/qrNbr read wallBoilingFlux_) and is credited
        // to this site in every lifecycle state, so the heater energy
        // balance does not depend on when the next bubble may form.
        scalar localBoilingPower = 0;
        scalar localConductionPower = 0;
        scalar localLiquidHeatingPower = 0;
        const bool siteActive =
            boilingClosure
         && (
                nearestSiteOverflow
              ? bool(prepassActive[siteI])
              : meanWallSuperheat + SMALL >= siteActivationSuperheat_[siteI]
            );
        // Below its activation superheat the site's area may still carry the
        // law, with all of its heat going to the liquid (inactiveAreaHeat).
        const bool inactiveToLiquid =
            boilingClosure && !siteActive && inactiveAreaHeat_ == "liquid";
        // liquidThenNearestSite(WithRoom): done before the loop; book it here.
        const bool inactiveFromPrepass = nearestSiteOverflow && !siteActive;
        const scalar evaporationFraction =
            inactiveToLiquid || inactiveFromPrepass
          ? scalar(0)
          : siteEvaporationFraction(meanWallSuperheat);
        if (inactiveFromPrepass)
        {
            // Only the liquid's share belongs to this site; what its liquid
            // could not take is booked at the active site that received it.
            localBoilingPower = inactiveLiquidPower[siteI];
            localLiquidHeatingPower = inactiveLiquidPower[siteI];
            localConductionPower = inactiveConductionPower[siteI];
            const scalar liquidEnergyStep = localLiquidHeatingPower*deltaT;
            cumulativeBoilingEnergy_[siteI] += liquidEnergyStep;
            cumulativeLiquidHeatingEnergy_[siteI] += liquidEnergyStep;
            totalBoilingPower += localBoilingPower;
            totalConductionPower += localConductionPower;
            totalLiquidHeatingPower += localLiquidHeatingPower;
        }
        else if (siteActive || inactiveToLiquid)
        {
            forAll(faces, faceListI)
            {
                const label faceI = faces[faceListI];
                const label celli = faceCells[faceI];
                const scalar area = mag(faceAreaVectors[faceI]);
                const scalar superheat =
                    max(wallTemperature[faceI] - TSat(celli), scalar(0));
                const scalar conductionFlux =
                    max(kappaPatch[faceI]*snGrad[faceI], scalar(0));
                scalar lawFlux =
                    wallBoilingCoefficient_
                   *Foam::pow(superheat, wallBoilingExponent_);
                if (wallBoilingMaxFlux_ > 0)
                {
                    lawFlux = min(lawFlux, wallBoilingMaxFlux_);
                }
                // "total": the law is the whole wall heat flux, so only
                // the part the mesh does not already conduct is added.
                scalar boilingFlux =
                    wallBoilingFluxMode_ == "total"
                  ? max(lawFlux - conductionFlux, scalar(0))
                  : lawFlux;
                // No nucleate boiling where the wall cell is dry.
                boilingFlux *=
                    max(min(phase1_[celli], scalar(1)), scalar(0));

                if (inactiveToLiquid)
                {
                    // Inactive area: the heat goes to the liquid of the wall
                    // cell, and only while that liquid is cooler than the
                    // wall. There is no bubble to send it to otherwise, so
                    // the face then carries no closure flux at all.
                    localConductionPower += conductionFlux*area;
                    if (liquidTemperature[celli] < wallTemperature[faceI])
                    {
                        const scalar liquidPower = boilingFlux*area;
                        boilingFluxPatch[faceI] = -boilingFlux;
                        wallLiquidHeat_[celli] += liquidPower/mesh.V()[celli];
                        localBoilingPower += liquidPower;
                        localLiquidHeatingPower += liquidPower;
                    }
                    continue;
                }

                boilingFluxPatch[faceI] = -boilingFlux;
                localBoilingPower += boilingFlux*area;
                localConductionPower += conductionFlux*area;

                // Wall heat partition: the share that does not form
                // bubbles at the wall heats the liquid of the wall cell,
                // but only while that liquid is cooler than the wall.
                if
                (
                    partitioned
                 && evaporationFraction < 1
                 && liquidTemperature[celli] < wallTemperature[faceI]
                )
                {
                    const scalar liquidPower =
                        (1 - evaporationFraction)*boilingFlux*area;
                    wallLiquidHeat_[celli] += liquidPower/mesh.V()[celli];
                    localLiquidHeatingPower += liquidPower;
                }
            }
            reduce(localBoilingPower, sumOp<scalar>());
            reduce(localConductionPower, sumOp<scalar>());
            if (partitioned || inactiveToLiquid)
            {
                reduce(localLiquidHeatingPower, sumOp<scalar>());
            }
            // Heat that the liquid over inactive sites could not take funds
            // this site's bubbles (liquidThenNearestSite(WithRoom)).
            if (nearestSiteOverflow)
            {
                localBoilingPower += receivedPower[siteI];
            }

            const scalar boilingEnergyStep = localBoilingPower*deltaT;
            const scalar liquidEnergyStep = localLiquidHeatingPower*deltaT;
            boilingEnergy_[siteI] += boilingEnergyStep - liquidEnergyStep;
            cumulativeBoilingEnergy_[siteI] += boilingEnergyStep;
            cumulativeLiquidHeatingEnergy_[siteI] += liquidEnergyStep;
            totalBoilingPower += localBoilingPower;
            totalConductionPower += localConductionPower;
            totalLiquidHeatingPower += localLiquidHeatingPower;
        }

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
        scalar localHandoffVaporVolume = 0;
        const labelList& releaseCells = sphereCells_[siteI];
        forAll(releaseCells, releaseI)
        {
            const label celli = releaseCells[releaseI];
            const scalar target = sphereTargetVapor_[siteI][releaseI];
            const scalar alphaLiquid =
                max(min(phase1_[celli], scalar(1)), scalar(0));
            const scalar alphaVapor = scalar(1) - alphaLiquid;
            localHandoffVaporVolume +=
                min(alphaVapor, target)*mesh.V()[celli];

            localSensibleEnergy +=
                sensibleReserveFraction_
               *phase1_[celli]
               *rhoL[celli]
               *cpL[celli]
               *mesh.V()[celli]
               *max(liquidTemperature[celli] - TSat(celli), scalar(0));
        }
        reduce(localSensibleEnergy, sumOp<scalar>());
        reduce(localHandoffVaporVolume, sumOp<scalar>());
        const scalar handoffVaporFraction =
            targetVaporVolume_[siteI] > SMALL
          ? min
            (
                localHandoffVaporVolume/targetVaporVolume_[siteI],
                scalar(1)
            )
          : scalar(0);
        const scalar macroCreatedVaporFraction =
            targetVaporVolume_[siteI] > SMALL
          ? min
            (
                createdVaporVolume_[siteI]/targetVaporVolume_[siteI],
                scalar(1)
            )
          : scalar(0);
        const scalar unfundedHandoffVaporFraction = max
        (
            handoffVaporFraction - macroCreatedVaporFraction,
            scalar(0)
        );
        const bool strictHandoffAccounting =
            maxPreexistingHandoffVaporFraction_ < scalar(1) - SMALL;

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
                initialSensibleCredit_[siteI] = localSensibleEnergy;
                chargedWallEnergy_[siteI] = 0.0;
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
                chargedWallEnergy_[siteI] += localWallEnergy;
                if
                (
                    storedEnergy_[siteI] + boilingEnergy_[siteI] + VSMALL
                 >= targetLatentEnergy_[siteI]
                )
                {
                    storedEnergy_[siteI] = min
                    (
                        storedEnergy_[siteI],
                        targetLatentEnergy_[siteI]
                    );

                    if
                    (
                        unfundedHandoffVaporFraction
                      <= maxPreexistingHandoffVaporFraction_ + SMALL
                    )
                    {
                        siteState_[siteI] = CREATING_DETACHED_BUBBLE;
                        creationStep_[siteI] = 0;
                        creationStartTime_[siteI] = timeValue;
                        if (writeHandoffAudit_ && Pstream::master())
                        {
                            Info<< "POPULATION_HANDOFF_FUNDED"
                                << " site=" << siteI
                                << " cycle=" << cycleId_[siteI]
                                << " time=" << timeValue
                                << " initialSensibleCredit="
                                << initialSensibleCredit_[siteI]
                                << " cumulativeWallEnergy="
                                << chargedWallEnergy_[siteI]
                                << " storedEnergy=" << storedEnergy_[siteI]
                                << " targetLatentEnergy="
                                << targetLatentEnergy_[siteI]
                                << " handoffVaporFraction="
                                << handoffVaporFraction
                                << endl;
                        }
                    }
                    else if
                    (
                        writeDiagnostics_
                     && Pstream::master()
                     &&
                        (
                            timeIndex % diagnosticsInterval_ == 0
                         || mesh.time().writeTime()
                        )
                    )
                    {
                        Info<< "POPULATION_HANDOFF_BLOCKED"
                            << " site=" << siteI
                            << " cycle=" << cycleId_[siteI]
                            << " time=" << timeValue
                            << " handoffVaporFraction="
                            << handoffVaporFraction
                            << " unfundedHandoffVaporFraction="
                            << unfundedHandoffVaporFraction
                            << " maximum="
                            << maxPreexistingHandoffVaporFraction_
                            << " storedEnergy=" << storedEnergy_[siteI]
                            << " targetEnergy="
                            << targetLatentEnergy_[siteI]
                            << endl;
                    }
                }
            }
        }
        else if (siteState_[siteI] == CREATING_DETACHED_BUBBLE)
        {
            const scalar handoffAge = timeValue - creationStartTime_[siteI];
            if
            (
                maximumHandoffDuration_ > 0
             && handoffAge > maximumHandoffDuration_ + SMALL
            )
            {
                FatalErrorInFunction
                    << "POPULATION_HANDOFF_TIMEOUT site=" << siteI
                    << " cycle=" << cycleId_[siteI]
                    << " time=" << timeValue
                    << " age=" << handoffAge
                    << " maximum=" << maximumHandoffDuration_
                    << " consumedLatentEnergy=" << consumedLatentEnergy_[siteI]
                    << " targetLatentEnergy=" << targetLatentEnergy_[siteI]
                    << " createdVaporVolume=" << createdVaporVolume_[siteI]
                    << " handoffVaporFraction=" << handoffVaporFraction
                    << " unfundedHandoffVaporFraction="
                    << unfundedHandoffVaporFraction
                    << exit(FatalError);
            }
            if (writeHandoffAudit_ && Pstream::master())
            {
                Info<< "POPULATION_HANDOFF_STATE"
                    << " site=" << siteI
                    << " cycle=" << cycleId_[siteI]
                    << " time=" << timeValue
                    << " age=" << handoffAge
                    << " handoffVaporFraction=" << handoffVaporFraction
                    << " macroCreatedVaporFraction=" << macroCreatedVaporFraction
                    << " unfundedHandoffVaporFraction="
                    << unfundedHandoffVaporFraction
                    << " storedEnergy=" << storedEnergy_[siteI]
                    << " consumedLatentEnergy=" << consumedLatentEnergy_[siteI]
                    << " targetLatentEnergy=" << targetLatentEnergy_[siteI]
                    << endl;
            }
            storedEnergy_[siteI] += localWallEnergy;
            if (strictHandoffAccounting)
            {
                storedEnergy_[siteI] = min
                (
                    storedEnergy_[siteI],
                    targetLatentEnergy_[siteI]
                );
            }
            creationStep_[siteI]++;

            scalar localRemainingVolume = 0;
            scalar localRequestedEnergy = 0;
            scalarField requestedAlpha(releaseCells.size(), 0.0);
            const label remainingSteps = max
            (
                creationSteps_ - creationStep_[siteI] + 1,
                label(1)
            );
            // Share of the still-missing vapor to request in this step.
            scalar stepFraction = scalar(1)/scalar(remainingSteps);
            if (creationDuration_ > 0)
            {
                const scalar remainingTime = creationDuration_ - handoffAge;
                stepFraction =
                    remainingTime > deltaT ? deltaT/remainingTime : scalar(1);
            }

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
                    remaining*stepFraction,
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
            const bool unfundedVaporBlocked =
                strictHandoffAccounting
             && unfundedHandoffVaporFraction
                > maxPreexistingHandoffVaporFraction_ + SMALL;
            const bool fullStencilBeforeBudget =
                strictHandoffAccounting
             && stencilComplete
             && !latentBudgetComplete;

            if (unfundedVaporBlocked || fullStencilBeforeBudget)
            {
                creationStep_[siteI] = max
                (
                    creationStep_[siteI] - 1,
                    label(0)
                );

                if
                (
                    writeDiagnostics_
                 && Pstream::master()
                 &&
                    (
                        timeIndex % diagnosticsInterval_ == 0
                     || mesh.time().writeTime()
                     || fullStencilBeforeBudget
                    )
                )
                {
                    Info<< "POPULATION_HANDOFF_BLOCKED"
                        << " site=" << siteI
                        << " cycle=" << cycleId_[siteI]
                        << " time=" << timeValue
                        << " stage=creating"
                        << " handoffVaporFraction="
                        << handoffVaporFraction
                        << " macroCreatedVaporFraction="
                        << macroCreatedVaporFraction
                        << " unfundedHandoffVaporFraction="
                        << unfundedHandoffVaporFraction
                        << " maximum="
                        << maxPreexistingHandoffVaporFraction_
                        << " stencilComplete=" << stencilComplete
                        << " remainingLatentBudget="
                        << remainingLatentBudget
                        << endl;
                }
            }
            else if
            (
                latentBudgetComplete
             || (!strictHandoffAccounting && stencilComplete)
            )
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
                        << " departureRadius=" << siteDepartureRadius_[siteI]
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
                        (storedEnergy_[siteI] + boilingEnergy_[siteI])
                       /localRequestedEnergy
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

                // Energy taken from the heater by the wall boiling closure
                // is spent first. That part of the bubble has already been
                // paid for by the solid, so it must not also be removed
                // from the fluid as a latent sink.
                const scalar requestedEnergy =
                    energyScale*localRequestedEnergy;
                const scalar boilingShare =
                    requestedEnergy > SMALL
                  ? min(boilingEnergy_[siteI], requestedEnergy)
                   /requestedEnergy
                  : scalar(0);

                scalar localCreatedVolume = 0;
                scalar localCreatedMass = 0;
                scalar localConsumedEnergy = 0;
                scalar localDisplacedVolumeRate = 0;
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
                    latentSink_[celli] +=
                        (scalar(1) - boilingShare)*latentRate;
                    activeReleaseMask_[celli] = scalar(1);
                    // Liquid in these cells is converted in place. Only
                    // rhoV/rhoL of it really evaporates; the rest has to
                    // be pushed aside, not deleted.
                    localDisplacedVolumeRate +=
                        alphaRate
                       *mesh.V()[celli]
                       *max(scalar(1) - rhoV[celli]/rhoL[celli], scalar(0));

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
                const scalar boilingConsumed =
                    boilingShare*localConsumedEnergy;
                boilingEnergy_[siteI] = max
                (
                    boilingEnergy_[siteI] - boilingConsumed,
                    scalar(0)
                );
                cumulativeBoilingBubbleEnergy_[siteI] += boilingConsumed;
                storedEnergy_[siteI] = max
                (
                    storedEnergy_[siteI]
                  - (localConsumedEnergy - boilingConsumed),
                    scalar(0)
                );
                cumulativeBirthVaporVolume_ += localCreatedVolume;

                if (continuitySourceMode_ == "displacedLiquid")
                {
                    // Return the displaced liquid as a volume source in
                    // the pure-liquid cells around the bubble. The solver
                    // divides the mass source by the mixture density, and
                    // its alpha equation keeps a full liquid cell full, so
                    // this adds exactly that volume of liquid.
                    reduce(localDisplacedVolumeRate, sumOp<scalar>());
                    const labelList& shellCells = displacementCells_[siteI];
                    const scalar liquidShellThreshold = 0.999;
                    scalar localShellVolume = 0;
                    forAll(shellCells, shellI)
                    {
                        const label celli = shellCells[shellI];
                        if (phase1_[celli] > liquidShellThreshold)
                        {
                            localShellVolume += mesh.V()[celli];
                        }
                    }
                    reduce(localShellVolume, sumOp<scalar>());

                    if (localShellVolume > SMALL)
                    {
                        const scalar volumeSource =
                            localDisplacedVolumeRate/localShellVolume;
                        forAll(shellCells, shellI)
                        {
                            const label celli = shellCells[shellI];
                            const scalar alphaLiquid = phase1_[celli];
                            if (alphaLiquid > liquidShellThreshold)
                            {
                                massBirthSource_[celli] +=
                                    volumeSource
                                   *(
                                        alphaLiquid*rhoL[celli]
                                      + (scalar(1) - alphaLiquid)*rhoV[celli]
                                    );
                            }
                        }
                        cumulativeDisplacedLiquidVolume_ +=
                            localDisplacedVolumeRate*deltaT;
                    }
                    else if
                    (
                        localDisplacedVolumeRate > SMALL
                     && Pstream::master()
                    )
                    {
                        Info<< "POPULATION_DISPLACEMENT_SKIPPED"
                            << " site=" << siteI
                            << " cycle=" << cycleId_[siteI]
                            << " time=" << timeValue
                            << " displacedVolumeRate="
                            << localDisplacedVolumeRate
                            << " reason=noLiquidShell"
                            << endl;
                    }
                }

                if (writeHandoffAudit_ && Pstream::master())
                {
                    Info<< "POPULATION_HANDOFF_SOURCE"
                        << " site=" << siteI
                        << " cycle=" << cycleId_[siteI]
                        << " time=" << timeValue
                        << " age=" << handoffAge
                        << " deltaVaporVolume=" << localCreatedVolume
                        << " deltaVaporMass=" << localCreatedMass
                        << " deltaLatentEnergy=" << localConsumedEnergy
                        << " createdVaporVolume=" << createdVaporVolume_[siteI]
                        << " createdVaporMass=" << createdVaporMass_[siteI]
                        << " consumedLatentEnergy=" << consumedLatentEnergy_[siteI]
                        << " targetLatentEnergy=" << targetLatentEnergy_[siteI]
                        << " remainingStencilVolume=" << localRemainingVolume
                        << endl;
                }

                const scalar postSourceRemainingLatentBudget = max
                (
                    targetLatentEnergy_[siteI]
                  - consumedLatentEnergy_[siteI],
                    scalar(0)
                );
                const bool creationOverdue =
                    creationDuration_ > 0
                  ? handoffAge > scalar(10)*creationDuration_
                  : creationStep_[siteI] >= maximumCreationSteps_;
                if
                (
                    creationOverdue
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
                >= rearmDelay(siteI)
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
                << " handoffVaporFraction=" << handoffVaporFraction
                << " unfundedHandoffVaporFraction="
                << unfundedHandoffVaporFraction
                << " captureArea=" << captureArea_[siteI]
                << " storedEnergy=" << storedEnergy_[siteI]
                << " targetEnergy=" << targetLatentEnergy_[siteI]
                << " createdVolume=" << createdVaporVolume_[siteI]
                << " targetVolume=" << targetVaporVolume_[siteI]
                << endl;

            if (boilingClosure)
            {
                Info<< "POPULATION_SITE_BOILING"
                    << " site=" << siteI
                    << " time=" << timeValue
                    << " meanWallSuperheat=" << meanWallSuperheat
                    << " boilingPower=" << localBoilingPower
                    << " conductionPower=" << localConductionPower
                    << " boilingEnergy=" << boilingEnergy_[siteI]
                    << " cumulativeBoilingEnergy="
                    << cumulativeBoilingEnergy_[siteI]
                    << " cumulativeBubbleEnergy="
                    << cumulativeBoilingBubbleEnergy_[siteI]
                    << " evaporationFraction=" << evaporationFraction
                    << " liquidHeatingPower=" << localLiquidHeatingPower
                    << " cumulativeLiquidHeatingEnergy="
                    << cumulativeLiquidHeatingEnergy_[siteI]
                    << " departureRadius=" << siteDepartureRadius_[siteI];
                if (nearestSiteOverflow)
                {
                    Info<< " receivedPower=" << receivedPower[siteI];
                }
                Info<< endl;
            }
        }
    }

    if (writeDiagnostics_ && diagnosticsStep && Pstream::master())
    {
        if (boilingClosure)
        {
            Info<< "POPULATION_BOILING_BALANCE"
                << " time=" << timeValue
                << " boilingPower=" << totalBoilingPower
                << " conductionPower=" << totalConductionPower
                << " cumulativeBoilingEnergy=" << sum(cumulativeBoilingEnergy_)
                << " cumulativeBubbleEnergy="
                << sum(cumulativeBoilingBubbleEnergy_)
                << " storedBoilingEnergy=" << sum(boilingEnergy_)
                << " liquidHeatingPower=" << totalLiquidHeatingPower
                << " cumulativeLiquidHeatingEnergy="
                << sum(cumulativeLiquidHeatingEnergy_);
            if (withRoom)
            {
                Info<< " overflowWithoutRoomPower=" << overflowWithoutRoomPower;
            }
            Info<< endl;
        }
        if (continuitySourceMode_ == "displacedLiquid")
        {
            Info<< "POPULATION_DISPLACEMENT_BALANCE"
                << " time=" << timeValue
                << " birthVaporVolume=" << cumulativeBirthVaporVolume_
                << " displacedLiquidVolume="
                << cumulativeDisplacedLiquidVolume_
                << endl;
        }
    }

    activeReleaseMask_.correctBoundaryConditions();
    siteStateField_.correctBoundaryConditions();
    alphaBirthSource_.correctBoundaryConditions();
    massBirthSource_.correctBoundaryConditions();
    latentSink_.correctBoundaryConditions();
    wallLiquidHeat_.correctBoundaryConditions();
}

void Foam::populationDetachedBubbleBirth::TSource1
(
    fvScalarMatrix& T1Eqn
)
{
    updateSources();
    // Latent heat of bubbles not paid for by the wall leaves the liquid; the
    // wall heat partition's share of the boiling heat enters it.
    T1Eqn.source() -=
        (latentSink_.internalField() - wallLiquidHeat_.internalField())
       *phase1_.mesh().V();
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