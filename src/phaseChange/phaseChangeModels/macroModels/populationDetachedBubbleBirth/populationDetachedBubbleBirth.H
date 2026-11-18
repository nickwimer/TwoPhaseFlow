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

    // Size of each released bubble. "fixed" (default) releases bubbles of
    // diameter 2*departureRadius. "logNormal" gives every release its own
    // diameter, drawn when the site rearms from a log-normal distribution
    // with median departureDiameterMedian and standard deviation of ln(D)
    // departureDiameterSigmaLn, redrawn until it lies in
    // [departureDiameterMin, departureDiameterMax]. The draw depends only on
    // (randomSeed, site, cycle), so every processor gets the same value
    // without communication. The site's release stencil and latent target
    // are rebuilt for that diameter, so the energy balance stays exact.
    word departureDiameterDistribution_;
    scalar departureDiameterMedian_;
    scalar departureDiameterSigmaLn_;
    scalar departureDiameterMin_;
    scalar departureDiameterMax_;
    label randomSeed_;

    // Release timing. rearmDelayScaling "diameter" scales minimumRearmDelay
    // by the last released diameter over 2*departureRadius, so a small
    // bubble clears its site sooner. waitingTimeMean > 0 adds an
    // exponentially distributed wait, drawn like the diameter, before a site
    // may activate again. Boiling energy keeps accumulating while a site
    // waits, so the scatter moves releases in time without losing energy.
    word rearmDelayScaling_;
    scalar waitingTimeMean_;

    // Share of an active site's boiling heat that funds its bubbles.
    // "vapor" (default): all of it, the original model. "constant": the
    // fraction evaporationFraction. "table": interpolated at the site's mean
    // wall superheat from evaporationFractionSuperheat/Value (clamped at the
    // ends). The rest heats the liquid in the wall cells of the site's
    // faces, but only while that liquid is cooler than the wall; otherwise
    // it too funds the site's bubbles, so the closure never heats the liquid
    // above the wall temperature and every joule still lands somewhere.
    word wallHeatPartition_;
    scalar evaporationFraction_;
    scalarField evaporationFractionSuperheat_;
    scalarField evaporationFractionValue_;

    // Heater area of sites below their activation superheat. "none"
    // (default): it carries no closure flux, the original model, so the wall
    // there warms until the site activates. With a fitted wall law (Cooper),
    // which already averages over the site density, that counts the density
    // twice. "liquid": the closure law applies there too and all of its heat
    // goes into the liquid of the wall cells, while that liquid is cooler
    // than the wall; a face whose liquid is not cooler carries no closure
    // flux. The inactive site's bubble energy is not charged.
    // "liquidThenNearestSite": the same, except that where the liquid is not
    // cooler the face keeps the law's flux and its heat funds the bubbles of
    // the nearest active site (in the heater plane), booked there in the same
    // step. So the whole heater always carries the law while any site is
    // active. These faces are handled in a pass before the site loop.
    word inactiveAreaHeat_;

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
    // Radius of the bubble each site will release next (or is releasing).
    scalarField siteDepartureRadius_;
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

    // Energy removed from the solid by the wall boiling closure that funds
    // bubbles and is not yet turned into vapor. It is never reset: every
    // joule of it must reappear as latent heat of a bubble.
    // cumulativeBoilingEnergy_ is all the closure has taken from the solid;
    // cumulativeLiquidHeatingEnergy_ is the part of it the wall heat
    // partition gave to the liquid instead, so that
    //     cumulativeBoilingEnergy = cumulativeBubbleEnergy
    //                             + boilingEnergy + cumulativeLiquidHeating.
    scalarField boilingEnergy_;
    scalarField cumulativeBoilingEnergy_;
    scalarField cumulativeBoilingBubbleEnergy_;
    scalarField cumulativeLiquidHeatingEnergy_;
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

    // Heat the wall heat partition puts into the liquid of the wall cells,
    // W/m3. Added to the liquid energy equation with the latent sink.
    volScalarField wallLiquidHeat_;

    enum siteState
    {
        WAITING_FOR_ACTIVATION = 0,
        CHARGING_DEPARTURE_BUBBLE = 1,
        CREATING_DETACHED_BUBBLE = 2,
        WAITING_FOR_CLEARANCE = 3
    };

    void validateControls() const;
    void buildSiteStencils();
    void buildReleaseStencil(const label siteI, const scalar radius);
    scalar hashUniform
    (
        const label siteI,
        const label cycle,
        const label stream
    ) const;
    scalar drawDepartureDiameter(const label siteI, const label cycle) const;
    scalar drawWaitingTime(const label siteI, const label cycle) const;
    scalar rearmDelay(const label siteI) const;
    scalar siteEvaporationFraction(const scalar superheat) const;
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
