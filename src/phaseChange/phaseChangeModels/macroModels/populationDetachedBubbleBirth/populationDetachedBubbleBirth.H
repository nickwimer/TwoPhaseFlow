/*---------------------------------------------------------------------------*\
    Coarse persistent-population boiling model with energy-limited detached-bubble release
\*---------------------------------------------------------------------------*/

#ifndef populationDetachedBubbleBirthMacro_H
#define populationDetachedBubbleBirthMacro_H

#include "typeInfo.H"
#include "macroModel.H"
#include "volFields.H"
#include "DynamicList.H"
#include "Switch.H"
#include "labelList.H"
#include "scalarField.H"

namespace Foam
{

class populationDetachedBubbleBirth
:
    public macroModel
{
private:

    word patch_;
    List<vector> siteCentres_;

    word energyOwnershipMode_;
    List<vector> energyOwnershipCentres_;
    labelList siteOwnershipIds_;

    scalar TSatValue_;
    scalarField siteActivationSuperheat_;
    scalar activationPersistence_;
    scalar activationStagger_;

    scalar captureRadius_;
    scalar departureRadius_;
    scalar detachmentGap_;
    scalar interfaceWidthCells_;
    scalar targetVaporFraction_;
    scalar maxAlphaVaporPerStep_;
    scalar maxPreexistingHandoffVaporFraction_;
    scalar maximumHandoffDuration_;
    label creationSteps_;
    label maximumCreationSteps_;

    // Time over which a funded bubble is created, in seconds. When positive
    // it replaces the step counts above, so the rate at which vapor is
    // created, and with displacedLiquid the speed of the liquid it pushes
    // aside, does not depend on the time step.
    scalar creationDuration_;

    scalar wallFluxUseFraction_;
    scalar sensibleReserveFraction_;

    scalar minimumRearmDelay_;
    scalar shieldingRadius_;
    scalar shieldingHeight_;
    scalar shieldingAlphaVaporThreshold_;

    word continuitySourceMode_;

    // Wall boiling closure. The resolved conduction flux into the first
    // liquid cell is limited by the wall-normal mesh spacing and is far below
    // a nucleate-boiling heat flux. With the closure enabled, each heater
    // face owned by an active site also passes the boiling heat flux
    //     qLaw = wallBoilingCoefficient*superheat^wallBoilingExponent
    // out of the solid. That flux is published as a boundary field which the
    // coupled temperature boundary conditions read through their qr/qrNbr
    // entries, and the same energy is credited to the owning site, where it
    // funds detached bubbles. Bubble formation itself is unchanged.
    word wallBoilingClosure_;
    scalar wallBoilingCoefficient_;
    scalar wallBoilingExponent_;
    scalar wallBoilingMaxFlux_;
    word wallBoilingFluxMode_;
    word wallBoilingFluxFieldName_;

    // Width, in local cell sizes, of the liquid shell around a new bubble
    // that receives the displaced liquid (continuitySourceMode
    // displacedLiquid).
    scalar displacementShellCells_;

    Switch writeDiagnostics_;
    Switch writeHandoffAudit_;
    word diagnosticPrefix_;
    label diagnosticsInterval_;

    label lastUpdateTimeIndex_;
    label patchID_;

    List<labelList> siteFaceLabels_;
    List<labelList> sphereCells_;
    List<scalarField> sphereTargetVapor_;
    List<labelList> shieldingCells_;
    List<labelList> displacementCells_;

    scalarField captureArea_;
    scalarField targetVaporVolume_;
    scalarField targetVaporMass_;
    scalarField targetLatentEnergy_;

    labelList siteState_;
    scalarField activationDuration_;
    scalarField storedEnergy_;
    scalarField releaseTime_;
    scalarField creationStartTime_;
    scalarField initialSensibleCredit_;
    scalarField chargedWallEnergy_;
    scalarField nextActivationTime_;
    labelList creationStep_;
    labelList cycleId_;

    scalarField createdVaporVolume_;
    scalarField createdVaporMass_;
    scalarField consumedLatentEnergy_;

    // Energy removed from the solid by the wall boiling closure and not yet
    // turned into vapor. It is never reset: every joule taken from the
    // heater must reappear as latent heat of a bubble.
    scalarField boilingEnergy_;
    scalarField cumulativeBoilingEnergy_;
    scalarField cumulativeBoilingBubbleEnergy_;
    scalar cumulativeDisplacedLiquidVolume_;
    scalar cumulativeBirthVaporVolume_;

    volScalarField activeReleaseMask_;
    volScalarField siteStateField_;
    volScalarField alphaBirthSource_;
    volScalarField massBirthSource_;
    volScalarField latentSink_;

    // Boundary value on the heater patch is the heat flux entering the
    // interface from the closure, in W/m2: minus the boiling flux.
    volScalarField wallBoilingFlux_;

    enum siteState
    {
        WAITING_FOR_ACTIVATION = 0,
        CHARGING_DEPARTURE_BUBBLE = 1,
        CREATING_DETACHED_BUBBLE = 2,
        WAITING_FOR_CLEARANCE = 3
    };

    void validateControls() const;
    void buildSiteStencils();
    scalar TSat(const label celli) const;
    const char* stateName(const label state) const;
    void resetSite(const label siteI, const scalar timeValue);
    void updateSources();

    populationDetachedBubbleBirth(const populationDetachedBubbleBirth&);
    void operator=(const populationDetachedBubbleBirth&);

public:

    TypeName("populationDetachedBubbleBirth");

    populationDetachedBubbleBirth
    (
        const phaseModel& phase1,
        const phaseModel& phase2,
        const volScalarField& p,
        singleComponentSatProp& satModel,
        const compressibleInterPhaseTransportModel& turbModel,
        const dictionary& dict
    );

    virtual ~populationDetachedBubbleBirth()
    {}

    virtual void momentumSource(fvVectorMatrix& UEqn);

    virtual void TSource1(fvScalarMatrix& T1Eqn);
    virtual void TSource2(fvScalarMatrix& T2Eqn);

    virtual void energySource(volScalarField& Q);
    virtual void energySource1(volScalarField& q1);
    virtual void energySource2(volScalarField& q2);

    virtual void massSource(volScalarField& rhoSource);
    virtual void alphaSource(volScalarField& rhoSource);
};

} // End namespace Foam

#endif

// ************************************************************************* //
