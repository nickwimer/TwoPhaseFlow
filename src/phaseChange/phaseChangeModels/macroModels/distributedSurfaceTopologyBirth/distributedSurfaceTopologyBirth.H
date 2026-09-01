/*---------------------------------------------------------------------------*\
    Distributed heated-surface topology-birth macro model
\*---------------------------------------------------------------------------*/

#ifndef distributedSurfaceTopologyBirthMacro_H
#define distributedSurfaceTopologyBirthMacro_H

#include "typeInfo.H"
#include "macroModel.H"
#include "volFields.H"
#include "DynamicList.H"
#include "Switch.H"

namespace Foam
{

class distributedSurfaceTopologyBirth
:
    public macroModel
{
private:

    wordList patches_;

    scalar TSatValue_;
    scalar superheatThreshold_;
    scalar baseHazardRate_;
    scalar hazardSuperheatScale_;
    scalar hazardExponent_;
    scalar maxHazardRate_;

    scalar minimumLiquidFraction_;

    scalar seedRadiusCells_;
    scalar maximumSeedRadius_;
    scalar seedCreationDuration_;
    scalar targetSeedVaporFraction_;
    scalar maxAlphaVaporPerStep_;
    scalar interfaceWidthCells_;
    scalar thermalReserveFraction_;
    word sourceCouplingMode_;

    scalar exclusionRadius_;
    scalar exclusionRadiusFactor_;
    scalar cooldownTime_;

    label maxBirthsPerStep_;
    label randomSeed_;

    Switch writeDiagnostics_;
    word diagnosticPrefix_;

    label lastUpdateTimeIndex_;
    label nextEventId_;

    scalar cumulativeTopologyMassDefect_;
    scalar cumulativeTopologyLatentDefect_;

    DynamicList<vector> eventCentres_;
    DynamicList<vector> eventNormals_;
    DynamicList<scalar> eventRadii_;
    DynamicList<scalar> eventEndTimes_;
    DynamicList<label> eventIds_;

    // Flat seed stencils are built once at event creation. This avoids
    // scanning the complete mesh for every active event on every CFD step.
    DynamicList<label> seedCellLabels_;
    DynamicList<scalar> seedCellTargetVapor_;
    DynamicList<scalar> seedCellEndTimes_;
    DynamicList<label> seedCellEventIds_;

    DynamicList<vector> recentCentres_;
    DynamicList<scalar> recentExclusionRadii_;
    DynamicList<scalar> recentExpiryTimes_;

    volScalarField superheat_;
    volScalarField eligibleMask_;
    volScalarField activeSeedMask_;
    volScalarField alphaBirthSource_;
    volScalarField massBirthSource_;
    volScalarField latentSink_;
    volScalarField equivalentMassDefectRate_;
    volScalarField equivalentLatentDefectRate_;

    void validateControls() const;

    scalar TSat(const label celli) const;

    scalar random01
    (
        const label timeIndex,
        const label patchOrdinal,
        const label facei
    ) const;

    bool recentlyExcluded
    (
        const vector& faceCentre,
        const scalar timeValue
    ) const;

    label activeEventCount(const scalar timeValue) const;

    label appendSeedStencil
    (
        const vector& wallCentre,
        const vector& inwardNormal,
        const scalar radius,
        const scalar eventEnd,
        const label eventId
    );

    label spawnEvents
    (
        const volScalarField& liquidTemperature
    );

    void applyActiveEvents
    (
        const volScalarField& liquidTemperature,
        const volScalarField& liquidDensity,
        const volScalarField& liquidCp,
        scalar& integratedAlphaSource,
        scalar& integratedMassSource,
        scalar& integratedLatentSink,
        scalar& integratedTopologyMassDefectRate,
        scalar& integratedTopologyLatentDefectRate,
        scalar& maxAlphaRate
    );

    void updateSources();

    distributedSurfaceTopologyBirth(const distributedSurfaceTopologyBirth&);
    void operator=(const distributedSurfaceTopologyBirth&);

public:

    TypeName("distributedSurfaceTopologyBirth");

    distributedSurfaceTopologyBirth
    (
        const phaseModel& phase1,
        const phaseModel& phase2,
        const volScalarField& p,
        singleComponentSatProp& satModel,
        const compressibleInterPhaseTransportModel& turbModel,
        const dictionary& dict
    );

    virtual ~distributedSurfaceTopologyBirth()
    {}

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
