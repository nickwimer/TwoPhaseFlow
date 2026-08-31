/*---------------------------------------------------------------------------*\
    Wall-superheat vapor-birth macro model for all-liquid heated-wall cases
\*---------------------------------------------------------------------------*/

#ifndef wallSuperheatVaporBirthMacro_H
#define wallSuperheatVaporBirthMacro_H

#include "typeInfo.H"
#include "macroModel.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "boolList.H"

namespace Foam
{

class wallSuperheatVaporBirth
:
    public macroModel
{
private:

    wordList patches_;
    scalar TSatValue_;
    scalar superheatOn_;
    scalar superheatOff_;
    scalar maxAlphaVaporPerStep_;
    label birthLayerCells_;
    scalar sourceRelaxation_;
    scalar epsLiquid_;
    Switch writeDiagnostics_;
    word diagnosticPrefix_;
    word sourceApplicationMode_;
    word postHandoffSourceMode_;
    word seedContinuitySourceMode_;
    word birthTriggerMode_;
    word seedHandoffMode_;
    word massDensityPhase_;
    word sourceFootprintMode_;
    word growthDepositionMode_;
    scalar growthBandHeight_;
    scalar growthBandAlphaMin_;
    scalar growthBandAlphaMax_;
    scalar growthBandNormalDecayLength_;
    scalar growthBandMinimumLiquid_;
    Switch growthBandRequirePrimaryBubble_;
    word siteMode_;
    word sitePatch_;
    vector siteCenter_;
    scalar siteRadius_;
    scalar contactLineAlphaMin_;
    scalar contactLineAlphaMax_;
    scalar minimumSourceLiquidFraction_;
    Switch contactLineNeighbourSearch_;
    scalar contactLineNeighbourDelta_;
    label maxActiveSites_;
    Switch lockSiteAfterBirth_;
    word sustainedBirthMode_;
    scalar minimumActiveDuration_;
    scalar minimumEvaporationSuperheat_;
    scalar thermalReserveFraction_;
    scalar minimumThermalAlphaStep_;
    word thermalLimitMode_;
    scalar wallFluxUseFraction_;
    scalar seedThermalReserveFraction_;
    scalar seedWallFluxUseFraction_;
    Switch writePerCellDiagnostics_;
    label perCellDiagnosticsInterval_;
    word perCellDiagnosticsFormat_;
    label perCellDiagnosticsSegmentID_;
    label lastPerCellDiagnosticsTimeIndex_;
    label compactDiagnosticsInterval_;
    mutable label lastStateDiagnosticsTimeIndex_;
    mutable label lastStateDiagnosticsState_;
    mutable label lastStateDiagnosticsCycleId_;
    mutable label lastSourceDiagnosticsTimeIndex_;
    mutable label lastSourceDiagnosticsState_;
    mutable label lastSourceDiagnosticsCycleId_;
    vector cavityCentre_;
    scalar cavityRadius_;
    scalar cavityContactAngle_;
    scalar cavitySurfaceTension_;
    label activationPersistenceSteps_;
    scalar activationPersistence_;
    label seedCreationSteps_;
    label maximumSeedCreationSteps_;
    scalar maxSeedAlphaVaporPerStep_;
    scalar targetSeedVaporFraction_;
    label minimumCellsPerHandoffRadius_;
    label preferredCellsPerHandoffRadius_;
    scalar maximumHandoffRadius_;
    Switch repeatedNucleationMode_;
    scalar departurePersistence_;
    scalar departureContactAreaThreshold_;
    scalar departureSiteVaporThreshold_;
    scalar rewetAlphaThreshold_;
    scalar rewetWallVaporThreshold_;
    scalar rewetPersistence_;
    scalar minimumRearmDelay_;
    Switch requireWettedSiteForActivation_;
    scalar activationAlphaLiquidThreshold_;
    scalar activationWallVaporThreshold_;
    scalar shieldingRadius_;
    scalar shieldingHeight_;
    scalar shieldingAlphaVaporThreshold_;
    scalar cavityCandidateRadius_;
    scalar cavityTipHeight_;
    scalar cavityCriticalTemperature_;
    scalar cavityCriticalSuperheat_;
    scalar handoffRadius_;
    scalar handoffCapHeight_;
    scalar handoffFootprintRadius_;
    scalar handoffTargetVolume_;
    scalar handoffTargetMass_;
    scalar handoffTargetLatentEnergy_;
    scalar localTangentialSpacing_;
    scalar localWallNormalSpacing_;

    labelList candidateCells_;
    labelList siteCells_;
    boolList siteCellMask_;
    scalarField wallFaceAreaCache_;
    labelList handoffCells_;
    boolList handoffCellMask_;
    scalarField handoffTargetVaporFractions_;
    labelList cavityTipCells_;
    scalarField cavityTipWeights_;
    labelList shieldingCells_;
    labelList sustainedState_;
    scalarField sustainedActivationTime_;
    label patchAdjacentCellCount_;
    label lastUpdateTimeIndex_;
    scalar firstTriggerTime_;
    scalar firstActiveSiteTime_;
    bool siteLocked_;
    label cavityBirthState_;
    label activationCounter_;
    scalar activationCandidateDuration_;
    label creationStepIndex_;
    scalar activationTime_;
    scalar handoffStartTime_;
    scalar handoffCompleteTime_;
    scalar storedHandoffEnergy_;
    scalar createdVaporVolume_;
    scalar createdVaporMass_;
    scalar displacedLiquidMass_;
    scalar consumedLatentEnergy_;
    scalar accumulatedWallEnergy_;
    scalar lastLocalSensibleEnergy_;
    scalar lastCavityTtip_;
    label cycleId_;
    label rearmCount_;
    scalar departureCandidateDuration_;
    scalar rewetCandidateDuration_;
    scalar departureTime_;
    scalar rewetTime_;
    scalar rearmTime_;
    scalar nextActivationTime_;
    bool departureDetected_;
    bool rewetDetected_;
    bool siteRearmed_;
    bool lastSiteShielded_;

    volScalarField superheat_;
    volScalarField activeMask_;
    volScalarField sustainedStateField_;
    volScalarField alphaBirthSource_;
    volScalarField massBirthSource_;
    volScalarField latentSink_;

    void validateControls() const;
    void buildCandidateCells();
    scalar siteVolume() const;
    scalar TSat(const label celli) const;
    scalar wallFluxIntoFluid(const label celli) const;
    scalar wallFaceArea(const label celli) const;
    bool contactLineEligible(const label celli) const;
    void appendInterfaceBandCells
    (
        const label wallCelli,
        const scalarField& alphaLiquid,
        boolList& included,
        scalarField& weights,
        scalarField& wallNormalDistance,
        label& eligibleCells
    ) const;
    void updateSources();
    bool cavityMode() const;
    const char* cavityStateName(const label state) const;
    void restoreCavityState();
    scalar interpolateCavityTipTemperature(const volScalarField& liquidTemperature) const;
    bool siteShielded() const;
    void buildHandoffCells(const volScalarField& vaporDensity);
    void updateCavityBirth
    (
        const volScalarField& liquidTemperature,
        const volScalarField& liquidDensity,
        const volScalarField& liquidCp,
        const volScalarField& vaporDensity
    );
    void writeCavityDiagnostics
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
        const bool writeSourceHistory = true
    ) const;
    //- Write only the source-history.csv row. Used from updateSources() during
    //  CONTACT_LINE_GROWTH when cavity diagnostics have already been written by
    //  updateCavityBirth() but source values are only known after the main loop.
    void writeSourceHistoryRow
    (
        const scalar integratedAlphaSource,
        const scalar integratedMassSource,
        const scalar integratedLatentSink,
        const label activeCells,
        const bool departureCandidate,
        const bool rewetCandidate
    ) const;
    void evaluateSiteAttachment
    (
        scalar& siteMeanAlphaLiquid,
        scalar& siteVaporWallFraction,
        scalar& attachedWallContactArea,
        bool& attachedBubblePresent
    ) const;
    void rearmForNextCycle();

    enum sustainedBirthState
    {
        INACTIVE = 0,
        REQUESTED_HOLD_ACTIVE = 1,
        THERMALLY_SUSPENDED = 2
    };

    enum cavityBirthState
    {
        WAITING_FOR_CRITERION = 0,
        ACCUMULATING_HANDOFF_ENERGY = 1,
        CREATING_RESOLVED_CAP = 2,
        CONTACT_LINE_GROWTH = 3,
        WAITING_FOR_REWET = 4
    };

    wallSuperheatVaporBirth(const wallSuperheatVaporBirth&);
    void operator=(const wallSuperheatVaporBirth&);

public:

    TypeName("wallSuperheatVaporBirth");

    wallSuperheatVaporBirth
    (
        const phaseModel& phase1,
        const phaseModel& phase2,
        const volScalarField& p,
        singleComponentSatProp& satModel,
        const compressibleInterPhaseTransportModel& turbModel,
        const dictionary& dict
    );

    virtual ~wallSuperheatVaporBirth()
    {}

    virtual void TSource1(fvScalarMatrix& T1Eqn);
    virtual void TSource2(fvScalarMatrix& T2Eqn);
    virtual void energySource(volScalarField& Q);
    virtual void energySource1(volScalarField& q1);
    virtual void energySource2(volScalarField& q2);
    virtual void massSource(volScalarField& rhoSource);
    virtual void alphaSource(volScalarField& rhoSource);

    static void resetSourceFields
    (
        volScalarField& superheat,
        volScalarField& alphaBirthSource,
        volScalarField& massBirthSource,
        volScalarField& latentSink
    );

    static bool dimensionSelfTest(const fvMesh& mesh, Ostream& os);

    static scalar thermallyBoundedAlphaStep
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
    );

    static scalar wallFluxAwareThermalAlphaStep
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
    );

    static labelList circularPatchSiteCells
    (
        const vectorField& faceCentres,
        const labelUList& faceCells,
        const vector& siteCenter,
        const scalar siteRadius
    );

    static bool contactLineEligible
    (
        const scalar alphaLiquid,
        const scalar neighbourMaxAlphaVapor,
        const scalar contactLineAlphaMin,
        const scalar contactLineAlphaMax,
        const scalar minimumSourceLiquidFraction,
        const bool contactLineNeighbourSearch,
        const scalar contactLineNeighbourDelta
    );
};

} // End namespace Foam

#endif

// ************************************************************************* //
