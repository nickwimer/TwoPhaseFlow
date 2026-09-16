/*---------------------------------------------------------------------------*\
    Coarse fixed-site boiling model with energy-limited detached-bubble release
\*---------------------------------------------------------------------------*/

#ifndef fixedSiteDetachedBubbleBirthMacro_H
#define fixedSiteDetachedBubbleBirthMacro_H

#include "typeInfo.H"
#include "macroModel.H"
#include "volFields.H"
#include "DynamicList.H"
#include "Switch.H"
#include "labelList.H"
#include "scalarField.H"

namespace Foam
{

class fixedSiteDetachedBubbleBirth
:
    public macroModel
{
private:

    word patch_;
    List<vector> siteCentres_;

    scalar TSatValue_;
    scalar activationSuperheat_;
    scalar activationPersistence_;
    scalar activationStagger_;

    scalar captureRadius_;
    scalar departureRadius_;
    scalar detachmentGap_;
    scalar interfaceWidthCells_;
    scalar targetVaporFraction_;
    scalar maxAlphaVaporPerStep_;
    label creationSteps_;
    label maximumCreationSteps_;

    scalar wallFluxUseFraction_;
    scalar sensibleReserveFraction_;

    scalar minimumRearmDelay_;
    scalar shieldingRadius_;
    scalar shieldingHeight_;
    scalar shieldingAlphaVaporThreshold_;

    word continuitySourceMode_;
    Switch writeDiagnostics_;
    word diagnosticPrefix_;
    label diagnosticsInterval_;

    label lastUpdateTimeIndex_;
    label patchID_;

    List<labelList> siteFaceLabels_;
    List<labelList> sphereCells_;
    List<scalarField> sphereTargetVapor_;
    List<labelList> shieldingCells_;

    scalarField captureArea_;
    scalarField targetVaporVolume_;
    scalarField targetVaporMass_;
    scalarField targetLatentEnergy_;

    labelList siteState_;
    scalarField activationDuration_;
    scalarField storedEnergy_;
    scalarField releaseTime_;
    scalarField nextActivationTime_;
    labelList creationStep_;
    labelList cycleId_;

    scalarField createdVaporVolume_;
    scalarField createdVaporMass_;
    scalarField consumedLatentEnergy_;

    volScalarField activeReleaseMask_;
    volScalarField siteStateField_;
    volScalarField alphaBirthSource_;
    volScalarField massBirthSource_;
    volScalarField latentSink_;

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

    fixedSiteDetachedBubbleBirth(const fixedSiteDetachedBubbleBirth&);
    void operator=(const fixedSiteDetachedBubbleBirth&);

public:

    TypeName("fixedSiteDetachedBubbleBirth");

    fixedSiteDetachedBubbleBirth
    (
        const phaseModel& phase1,
        const phaseModel& phase2,
        const volScalarField& p,
        singleComponentSatProp& satModel,
        const compressibleInterPhaseTransportModel& turbModel,
        const dictionary& dict
    );

    virtual ~fixedSiteDetachedBubbleBirth()
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
