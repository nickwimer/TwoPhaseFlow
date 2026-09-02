/*---------------------------------------------------------------------------*\
    Distributed heated-surface topology-birth macro model
\*---------------------------------------------------------------------------*/

#include "distributedSurfaceTopologyBirth.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvPatch.H"
#include "fvScalarMatrix.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"

#include <cmath>
#include <cstdint>

namespace Foam
{

defineTypeNameAndDebug(distributedSurfaceTopologyBirth, 0);
addToRunTimeSelectionTable
(
    macroModel,
    distributedSurfaceTopologyBirth,
    components
);


distributedSurfaceTopologyBirth::distributedSurfaceTopologyBirth
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
    superheatThreshold_(modelDict().get<scalar>("superheatThreshold")),
    baseHazardRate_(modelDict().get<scalar>("baseHazardRate")),
    hazardSuperheatScale_
    (
        modelDict().lookupOrDefault<scalar>("hazardSuperheatScale", 1.0)
    ),
    hazardExponent_
    (
        modelDict().lookupOrDefault<scalar>("hazardExponent", 1.0)
    ),
    maxHazardRate_
    (
        modelDict().lookupOrDefault<scalar>("maxHazardRate", 0.0)
    ),
    minimumLiquidFraction_
    (
        modelDict().lookupOrDefault<scalar>("minimumLiquidFraction", 0.98)
    ),
    seedRadiusCells_
    (
        modelDict().lookupOrDefault<scalar>("seedRadiusCells", 3.0)
    ),
    maximumSeedRadius_
    (
        modelDict().lookupOrDefault<scalar>("maximumSeedRadius", 0.0)
    ),
    seedCreationDuration_
    (
        modelDict().lookupOrDefault<scalar>("seedCreationDuration", 1.0e-4)
    ),
    targetSeedVaporFraction_
    (
        modelDict().lookupOrDefault<scalar>("targetSeedVaporFraction", 0.9)
    ),
    maxAlphaVaporPerStep_
    (
        modelDict().lookupOrDefault<scalar>("maxAlphaVaporPerStep", 0.1)
    ),
    interfaceWidthCells_
    (
        modelDict().lookupOrDefault<scalar>("interfaceWidthCells", 1.0)
    ),
    thermalReserveFraction_
    (
        modelDict().lookupOrDefault<scalar>("thermalReserveFraction", 0.0)
    ),
    sourceCouplingMode_
    (
        modelDict().lookupOrDefault<word>
        (
            "sourceCouplingMode",
            "conservative"
        )
    ),
    exclusionRadius_
    (
        modelDict().lookupOrDefault<scalar>("exclusionRadius", 0.0)
    ),
    exclusionRadiusFactor_
    (
        modelDict().lookupOrDefault<scalar>("exclusionRadiusFactor", 2.5)
    ),
    cooldownTime_
    (
        modelDict().lookupOrDefault<scalar>("cooldownTime", 0.0)
    ),
    maxBirthsPerStep_
    (
        modelDict().lookupOrDefault<label>("maxBirthsPerStep", 4)
    ),
    randomSeed_
    (
        modelDict().lookupOrDefault<label>("randomSeed", 20401)
    ),
    writeDiagnostics_
    (
        modelDict().lookupOrDefault<Switch>("writeDiagnostics", true)
    ),
    diagnosticPrefix_
    (
        modelDict().lookupOrDefault<word>("diagnosticPrefix", typeName)
    ),
    lastUpdateTimeIndex_(-1),
    nextEventId_(0),
    cumulativeTopologyMassDefect_(0),
    cumulativeTopologyLatentDefect_(0),
    eventCentres_(),
    eventNormals_(),
    eventRadii_(),
    eventEndTimes_(),
    eventIds_(),
    seedCellLabels_(),
    seedCellTargetVapor_(),
    seedCellEndTimes_(),
    seedCellEventIds_(),
    recentCentres_(),
    recentExclusionRadii_(),
    recentExpiryTimes_(),
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
        dimensionedScalar("zero", dimTemperature, 0.0),
        "zeroGradient"
    ),
    eligibleMask_
    (
        IOobject
        (
            word(diagnosticPrefix_ + "EligibleMask"),
            phase1.mesh().time().timeName(),
            phase1.mesh(),
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        phase1.mesh(),
        dimensionedScalar("zero", dimless, 0.0),
        "zeroGradient"
    ),
    activeSeedMask_
    (
        IOobject
        (
            word(diagnosticPrefix_ + "ActiveSeedMask"),
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
    equivalentMassDefectRate_
    (
        IOobject
        (
            word(diagnosticPrefix_ + "EquivalentMassDefectRate"),
            phase1.mesh().time().timeName(),
            phase1.mesh(),
            IOobject::NO_READ,
            writeDiagnostics_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        phase1.mesh(),
        dimensionedScalar("zero", dimDensity/dimTime, 0.0),
        "zeroGradient"
    ),
    equivalentLatentDefectRate_
    (
        IOobject
        (
            word(diagnosticPrefix_ + "EquivalentLatentDefectRate"),
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
    validateControls();

    if (Pstream::master())
    {
        Info<< "distributedSurfaceTopologyBirth: patches=" << patches_
            << ", superheatThreshold=" << superheatThreshold_
            << " K, baseHazardRate=" << baseHazardRate_
            << " 1/(m2 s), seedRadiusCells=" << seedRadiusCells_
            << ", seedCreationDuration=" << seedCreationDuration_ << " s"
            << ", targetSeedVaporFraction=" << targetSeedVaporFraction_
            << ", sourceCouplingMode=" << sourceCouplingMode_
            << ", parallel=" << Pstream::parRun()
            << endl;
    }
}


void distributedSurfaceTopologyBirth::validateControls() const
{
    if (patches_.empty())
    {
        FatalErrorInFunction
            << "At least one heated-wall patch is required."
            << exit(FatalError);
    }

    const fvMesh& mesh = phase1_.mesh();
    forAll(patches_, patchi)
    {
        if (mesh.boundaryMesh().findPatchID(patches_[patchi]) < 0)
        {
            FatalErrorInFunction
                << "Cannot find patch " << patches_[patchi]
                << " in region " << mesh.name()
                << exit(FatalError);
        }
    }

    if (baseHazardRate_ < 0)
    {
        FatalErrorInFunction
            << "baseHazardRate must be non-negative."
            << exit(FatalError);
    }

    if (hazardSuperheatScale_ <= SMALL)
    {
        FatalErrorInFunction
            << "hazardSuperheatScale must be positive."
            << exit(FatalError);
    }

    if (hazardExponent_ < 0)
    {
        FatalErrorInFunction
            << "hazardExponent must be non-negative."
            << exit(FatalError);
    }

    if (maxHazardRate_ < 0)
    {
        FatalErrorInFunction
            << "maxHazardRate must be non-negative; use 0 for no cap."
            << exit(FatalError);
    }

    if
    (
        minimumLiquidFraction_ < 0
     || minimumLiquidFraction_ > 1
    )
    {
        FatalErrorInFunction
            << "minimumLiquidFraction must be in [0,1]."
            << exit(FatalError);
    }

    if (seedRadiusCells_ <= SMALL)
    {
        FatalErrorInFunction
            << "seedRadiusCells must be positive."
            << exit(FatalError);
    }

    if (maximumSeedRadius_ < 0)
    {
        FatalErrorInFunction
            << "maximumSeedRadius must be non-negative; use 0 for no cap."
            << exit(FatalError);
    }

    if (seedCreationDuration_ <= SMALL)
    {
        FatalErrorInFunction
            << "seedCreationDuration must be positive."
            << exit(FatalError);
    }

    if
    (
        targetSeedVaporFraction_ <= 0
     || targetSeedVaporFraction_ >= 1
    )
    {
        FatalErrorInFunction
            << "targetSeedVaporFraction must be in (0,1)."
            << exit(FatalError);
    }

    if
    (
        maxAlphaVaporPerStep_ <= 0
     || maxAlphaVaporPerStep_ > 1
    )
    {
        FatalErrorInFunction
            << "maxAlphaVaporPerStep must be in (0,1]."
            << exit(FatalError);
    }

    if (interfaceWidthCells_ <= SMALL)
    {
        FatalErrorInFunction
            << "interfaceWidthCells must be positive."
            << exit(FatalError);
    }

    if
    (
        thermalReserveFraction_ < 0
     || thermalReserveFraction_ > 1
    )
    {
        FatalErrorInFunction
            << "thermalReserveFraction must be in [0,1]."
            << exit(FatalError);
    }

    if
    (
        sourceCouplingMode_ != "conservative"
     && sourceCouplingMode_ != "topologyOnly"
    )
    {
        FatalErrorInFunction
            << "sourceCouplingMode must be conservative or topologyOnly; found "
            << sourceCouplingMode_
            << exit(FatalError);
    }

    if
    (
        exclusionRadius_ < 0
     || exclusionRadiusFactor_ <= SMALL
     || cooldownTime_ < 0
    )
    {
        FatalErrorInFunction
            << "exclusionRadius and cooldownTime must be non-negative and "
            << "exclusionRadiusFactor must be positive."
            << exit(FatalError);
    }

    if (maxBirthsPerStep_ < 1)
    {
        FatalErrorInFunction
            << "maxBirthsPerStep must be >= 1."
            << exit(FatalError);
    }
}


scalar distributedSurfaceTopologyBirth::TSat(const label celli) const
{
    return TSatValue_ > 0 ? TSatValue_ : satModel_.TSat()[celli];
}


scalar distributedSurfaceTopologyBirth::random01
(
    const label timeIndex,
    const label patchOrdinal,
    const vector& faceCentre
) const
{
    // Quantize geometry to 1 nm before hashing. This is many orders of
    // magnitude finer than the capability mesh, while avoiding decomposition-
    // dependent processor-local face numbering and tiny coordinate roundoff.
    const std::int64_t qx =
        static_cast<std::int64_t>(std::llround(faceCentre.x()*1.0e9));
    const std::int64_t qy =
        static_cast<std::int64_t>(std::llround(faceCentre.y()*1.0e9));
    const std::int64_t qz =
        static_cast<std::int64_t>(std::llround(faceCentre.z()*1.0e9));

    // Zig-zag signed coordinates before combining them. Directly converting
    // two's-complement signed coordinates to uint64_t and XORing the terms
    // gives exact inversion symmetry for opposite heater locations.
    const auto signedCoordinateKey = [](const std::int64_t q)
    {
        const std::uint64_t magnitude =
            q >= 0
          ? static_cast<std::uint64_t>(q)
          : static_cast<std::uint64_t>(-(q + 1)) + UINT64_C(1);

        return
            q >= 0
          ? (magnitude << 1)
          : ((magnitude << 1) - UINT64_C(1));
    };

    const std::uint64_t kx = signedCoordinateKey(qx);
    const std::uint64_t ky = signedCoordinateKey(qy);
    const std::uint64_t kz = signedCoordinateKey(qz);

    std::uint64_t x =
        static_cast<std::uint64_t>(static_cast<std::int64_t>(randomSeed_));

    x ^= static_cast<std::uint64_t>(timeIndex + 1)
       * UINT64_C(0x9E3779B97F4A7C15);
    x ^= static_cast<std::uint64_t>(patchOrdinal + 1)
       * UINT64_C(0xBF58476D1CE4E5B9);
    x ^= kx * UINT64_C(0x94D049BB133111EB);
    x ^= ky * UINT64_C(0xD6E8FEB86659FD93);
    x ^= kz * UINT64_C(0xA5A3564E27F8862D);

    x ^= x >> 30;
    x *= UINT64_C(0xBF58476D1CE4E5B9);
    x ^= x >> 27;
    x *= UINT64_C(0x94D049BB133111EB);
    x ^= x >> 31;

    const std::uint64_t mantissa = x >> 11;
    return scalar(mantissa)*(scalar(1.0)/scalar(9007199254740992.0));
}


bool distributedSurfaceTopologyBirth::recentlyExcluded
(
    const vector& faceCentre,
    const scalar timeValue
) const
{
    forAll(recentCentres_, i)
    {
        if
        (
            recentExpiryTimes_[i] + SMALL >= timeValue
         && magSqr(faceCentre - recentCentres_[i])
            < sqr(recentExclusionRadii_[i])
        )
        {
            return true;
        }
    }

    return false;
}


label distributedSurfaceTopologyBirth::activeEventCount
(
    const scalar timeValue
) const
{
    label count = 0;

    forAll(eventEndTimes_, i)
    {
        if (eventEndTimes_[i] + SMALL >= timeValue)
        {
            ++count;
        }
    }

    return count;
}


label distributedSurfaceTopologyBirth::appendSeedStencil
(
    const vector& wallCentre,
    const vector& inwardNormal,
    const scalar radius,
    const scalar eventEnd,
    const label eventId
)
{
    const fvMesh& mesh = phase1_.mesh();
    const volVectorField& cellCentres = mesh.C();

    label appended = 0;

    forAll(phase1_, celli)
    {
        const scalar cellVolume = mesh.V()[celli];
        if (cellVolume <= VSMALL)
        {
            continue;
        }

        const scalar cellSize =
            std::cbrt(max(cellVolume, scalar(VSMALL)));
        const vector displacement =
            cellCentres[celli] - wallCentre;
        const scalar wallNormalDistance =
            displacement & inwardNormal;

        if
        (
            wallNormalDistance < -0.5*cellSize
         || wallNormalDistance > radius + cellSize
        )
        {
            continue;
        }

        const scalar interfaceWidth =
            max(interfaceWidthCells_*cellSize, scalar(SMALL));
        const scalar radialDistance = mag(displacement);

        scalar shape =
            (radius + 0.5*interfaceWidth - radialDistance)
           /interfaceWidth;
        shape = min(max(shape, scalar(0)), scalar(1));

        if (shape <= SMALL)
        {
            continue;
        }

        seedCellLabels_.append(celli);
        seedCellTargetVapor_.append(targetSeedVaporFraction_*shape);
        seedCellEndTimes_.append(eventEnd);
        seedCellEventIds_.append(eventId);
        ++appended;
    }

    return appended;
}


label distributedSurfaceTopologyBirth::spawnEvents
(
    const volScalarField& liquidTemperature
)
{
    const fvMesh& mesh = phase1_.mesh();
    const scalar deltaT = mesh.time().deltaTValue();
    const scalar timeValue = mesh.time().value();
    const label timeIndex = mesh.time().timeIndex();

    if (deltaT <= SMALL || baseHazardRate_ <= 0)
    {
        return 0;
    }

    DynamicList<vector> candidateCentres;
    DynamicList<vector> candidateNormals;
    DynamicList<scalar> candidateRadii;
    DynamicList<scalar> candidateSuperheats;
    DynamicList<scalar> candidateHazards;
    DynamicList<scalar> candidateProbabilities;
    DynamicList<scalar> candidateExclusionRadii;
    DynamicList<label> candidatePatchOrdinals;
    DynamicList<label> candidateLocalFaces;

    forAll(patches_, patchOrdinal)
    {
        const label patchID =
            mesh.boundaryMesh().findPatchID(patches_[patchOrdinal]);
        const fvPatch& patch = mesh.boundary()[patchID];
        const labelUList& faceCells = patch.faceCells();
        const vectorField& faceCentres = patch.Cf();
        const vectorField& faceAreaVectors = patch.Sf();
        const scalarField& wallTemperature =
            liquidTemperature.boundaryField()[patchID];

        forAll(faceCells, facei)
        {
            const label celli = faceCells[facei];
            const scalar wallSuperheat =
                wallTemperature[facei] - TSat(celli);

            superheat_[celli] =
                max(superheat_[celli], wallSuperheat);

            const bool thermallyEligible =
                wallSuperheat > superheatThreshold_;
            const bool liquidEligible =
                phase1_[celli] >= minimumLiquidFraction_;

            if
            (
                !thermallyEligible
             || !liquidEligible
             || recentlyExcluded(faceCentres[facei], timeValue)
            )
            {
                continue;
            }

            eligibleMask_[celli] = 1.0;

            const scalar normalizedExcess =
                max
                (
                    (wallSuperheat - superheatThreshold_)
                   /hazardSuperheatScale_,
                    scalar(0)
                );

            scalar hazard =
                baseHazardRate_
               *std::pow(normalizedExcess, hazardExponent_);

            if (maxHazardRate_ > 0)
            {
                hazard = min(hazard, maxHazardRate_);
            }

            const scalar faceArea = mag(faceAreaVectors[facei]);
            const scalar eventProbability =
                min
                (
                    max
                    (
                        scalar(1)
                      - std::exp(-hazard*faceArea*deltaT),
                        scalar(0)
                    ),
                    scalar(1)
                );

            if
            (
                random01(timeIndex, patchOrdinal, faceCentres[facei])
              >= eventProbability
            )
            {
                continue;
            }

            const scalar localDelta =
                std::cbrt(max(mesh.V()[celli], scalar(VSMALL)));

            scalar radius = seedRadiusCells_*localDelta;
            if (maximumSeedRadius_ > 0)
            {
                radius = min(radius, maximumSeedRadius_);
            }

            const scalar faceAreaMag =
                max(faceArea, scalar(VSMALL));
            const vector inwardNormal =
                -faceAreaVectors[facei]/faceAreaMag;
            const scalar localExclusionRadius =
                exclusionRadius_ > 0
              ? exclusionRadius_
              : exclusionRadiusFactor_*radius;

            candidateCentres.append(faceCentres[facei]);
            candidateNormals.append(inwardNormal);
            candidateRadii.append(radius);
            candidateSuperheats.append(wallSuperheat);
            candidateHazards.append(hazard);
            candidateProbabilities.append(eventProbability);
            candidateExclusionRadii.append(localExclusionRadius);
            candidatePatchOrdinals.append(patchOrdinal);
            candidateLocalFaces.append(facei);
        }
    }

    const List<vector> localCentres(candidateCentres);
    const List<vector> localNormals(candidateNormals);
    const List<scalar> localRadii(candidateRadii);
    const List<scalar> localSuperheats(candidateSuperheats);
    const List<scalar> localHazards(candidateHazards);
    const List<scalar> localProbabilities(candidateProbabilities);
    const List<scalar> localExclusionRadii(candidateExclusionRadii);
    const List<label> localPatchOrdinals(candidatePatchOrdinals);
    const List<label> localFaces(candidateLocalFaces);

    const List<List<vector>> allCentres =
        Pstream::listGatherValues(localCentres);
    const List<List<vector>> allNormals =
        Pstream::listGatherValues(localNormals);
    const List<List<scalar>> allRadii =
        Pstream::listGatherValues(localRadii);
    const List<List<scalar>> allSuperheats =
        Pstream::listGatherValues(localSuperheats);
    const List<List<scalar>> allHazards =
        Pstream::listGatherValues(localHazards);
    const List<List<scalar>> allProbabilities =
        Pstream::listGatherValues(localProbabilities);
    const List<List<scalar>> allExclusionRadii =
        Pstream::listGatherValues(localExclusionRadii);
    const List<List<label>> allPatchOrdinals =
        Pstream::listGatherValues(localPatchOrdinals);
    const List<List<label>> allFaces =
        Pstream::listGatherValues(localFaces);

    List<vector> acceptedCentres;
    List<vector> acceptedNormals;
    List<scalar> acceptedRadii;
    List<scalar> acceptedSuperheats;
    List<scalar> acceptedHazards;
    List<scalar> acceptedProbabilities;
    List<scalar> acceptedExclusionRadii;
    List<label> acceptedPatchOrdinals;
    List<label> acceptedSourceProcs;
    List<label> acceptedSourceFaces;
    List<label> acceptedEventIds;

    if (Pstream::master())
    {
        label totalCandidates = 0;
        forAll(allCentres, proci)
        {
            totalCandidates += allCentres[proci].size();
        }

        List<vector> flatCentres(totalCandidates);
        List<vector> flatNormals(totalCandidates);
        List<scalar> flatRadii(totalCandidates);
        List<scalar> flatSuperheats(totalCandidates);
        List<scalar> flatHazards(totalCandidates);
        List<scalar> flatProbabilities(totalCandidates);
        List<scalar> flatExclusionRadii(totalCandidates);
        List<label> flatPatchOrdinals(totalCandidates);
        List<label> flatSourceProcs(totalCandidates);
        List<label> flatSourceFaces(totalCandidates);
        List<label> order(totalCandidates);

        label flatI = 0;
        forAll(allCentres, proci)
        {
            forAll(allCentres[proci], locali)
            {
                flatCentres[flatI] = allCentres[proci][locali];
                flatNormals[flatI] = allNormals[proci][locali];
                flatRadii[flatI] = allRadii[proci][locali];
                flatSuperheats[flatI] = allSuperheats[proci][locali];
                flatHazards[flatI] = allHazards[proci][locali];
                flatProbabilities[flatI] = allProbabilities[proci][locali];
                flatExclusionRadii[flatI] =
                    allExclusionRadii[proci][locali];
                flatPatchOrdinals[flatI] =
                    allPatchOrdinals[proci][locali];
                flatSourceProcs[flatI] = proci;
                flatSourceFaces[flatI] = allFaces[proci][locali];
                order[flatI] = flatI;
                ++flatI;
            }
        }

        // Stable physical ordering, independent of processor decomposition.
        // The face-centre coordinates are unique on a given configured patch.
        for (label i = 1; i < order.size(); ++i)
        {
            const label key = order[i];
            label j = i - 1;

            auto before = [&](const label a, const label b)
            {
                if (flatPatchOrdinals[a] != flatPatchOrdinals[b])
                {
                    return flatPatchOrdinals[a] < flatPatchOrdinals[b];
                }

                const vector& ca = flatCentres[a];
                const vector& cb = flatCentres[b];

                if (ca.x() != cb.x()) return ca.x() < cb.x();
                if (ca.y() != cb.y()) return ca.y() < cb.y();
                if (ca.z() != cb.z()) return ca.z() < cb.z();
                if (flatSourceProcs[a] != flatSourceProcs[b])
                {
                    return flatSourceProcs[a] < flatSourceProcs[b];
                }
                return flatSourceFaces[a] < flatSourceFaces[b];
            };

            while (j >= 0 && before(key, order[j]))
            {
                order[j + 1] = order[j];
                --j;
            }
            order[j + 1] = key;
        }

        acceptedCentres.resize(maxBirthsPerStep_);
        acceptedNormals.resize(maxBirthsPerStep_);
        acceptedRadii.resize(maxBirthsPerStep_);
        acceptedSuperheats.resize(maxBirthsPerStep_);
        acceptedHazards.resize(maxBirthsPerStep_);
        acceptedProbabilities.resize(maxBirthsPerStep_);
        acceptedExclusionRadii.resize(maxBirthsPerStep_);
        acceptedPatchOrdinals.resize(maxBirthsPerStep_);
        acceptedSourceProcs.resize(maxBirthsPerStep_);
        acceptedSourceFaces.resize(maxBirthsPerStep_);
        acceptedEventIds.resize(maxBirthsPerStep_);

        label acceptedCount = 0;

        forAll(order, orderI)
        {
            if (acceptedCount >= maxBirthsPerStep_)
            {
                break;
            }

            const label candidateI = order[orderI];
            const vector& centre = flatCentres[candidateI];

            bool excluded = recentlyExcluded(centre, timeValue);

            for
            (
                label acceptedI = 0;
                !excluded && acceptedI < acceptedCount;
                ++acceptedI
            )
            {
                if
                (
                    magSqr(centre - acceptedCentres[acceptedI])
                  < sqr(acceptedExclusionRadii[acceptedI])
                )
                {
                    excluded = true;
                }
            }

            if (excluded)
            {
                continue;
            }

            acceptedCentres[acceptedCount] = centre;
            acceptedNormals[acceptedCount] = flatNormals[candidateI];
            acceptedRadii[acceptedCount] = flatRadii[candidateI];
            acceptedSuperheats[acceptedCount] = flatSuperheats[candidateI];
            acceptedHazards[acceptedCount] = flatHazards[candidateI];
            acceptedProbabilities[acceptedCount] =
                flatProbabilities[candidateI];
            acceptedExclusionRadii[acceptedCount] =
                flatExclusionRadii[candidateI];
            acceptedPatchOrdinals[acceptedCount] =
                flatPatchOrdinals[candidateI];
            acceptedSourceProcs[acceptedCount] =
                flatSourceProcs[candidateI];
            acceptedSourceFaces[acceptedCount] =
                flatSourceFaces[candidateI];
            acceptedEventIds[acceptedCount] = nextEventId_++;
            ++acceptedCount;
        }

        acceptedCentres.resize(acceptedCount);
        acceptedNormals.resize(acceptedCount);
        acceptedRadii.resize(acceptedCount);
        acceptedSuperheats.resize(acceptedCount);
        acceptedHazards.resize(acceptedCount);
        acceptedProbabilities.resize(acceptedCount);
        acceptedExclusionRadii.resize(acceptedCount);
        acceptedPatchOrdinals.resize(acceptedCount);
        acceptedSourceProcs.resize(acceptedCount);
        acceptedSourceFaces.resize(acceptedCount);
        acceptedEventIds.resize(acceptedCount);
    }

    Pstream::broadcastList(acceptedCentres);
    Pstream::broadcastList(acceptedNormals);
    Pstream::broadcastList(acceptedRadii);
    Pstream::broadcastList(acceptedSuperheats);
    Pstream::broadcastList(acceptedHazards);
    Pstream::broadcastList(acceptedProbabilities);
    Pstream::broadcastList(acceptedExclusionRadii);
    Pstream::broadcastList(acceptedPatchOrdinals);
    Pstream::broadcastList(acceptedSourceProcs);
    Pstream::broadcastList(acceptedSourceFaces);
    Pstream::broadcastList(acceptedEventIds);

    const scalar eventEnd = timeValue + seedCreationDuration_;

    forAll(acceptedCentres, eventI)
    {
        const label eventId = acceptedEventIds[eventI];

        if (!Pstream::master())
        {
            nextEventId_ = max(nextEventId_, eventId + 1);
        }

        eventCentres_.append(acceptedCentres[eventI]);
        eventNormals_.append(acceptedNormals[eventI]);
        eventRadii_.append(acceptedRadii[eventI]);
        eventEndTimes_.append(eventEnd);
        eventIds_.append(eventId);

        label stencilCells =
            appendSeedStencil
            (
                acceptedCentres[eventI],
                acceptedNormals[eventI],
                acceptedRadii[eventI],
                eventEnd,
                eventId
            );

        reduce(stencilCells, sumOp<label>());

        if (stencilCells == 0)
        {
            FatalErrorInFunction
                << "Event " << eventId
                << " generated an empty global seed stencil at "
                << acceptedCentres[eventI]
                << exit(FatalError);
        }

        recentCentres_.append(acceptedCentres[eventI]);
        recentExclusionRadii_.append(acceptedExclusionRadii[eventI]);
        recentExpiryTimes_.append(eventEnd + cooldownTime_);

        if (Pstream::master())
        {
            const label patchOrdinal = acceptedPatchOrdinals[eventI];

            Info<< "DSTB_EVENT"
                << " eventId=" << eventId
                << " time=" << timeValue
                << " patch=" << patches_[patchOrdinal]
                << " face=" << acceptedSourceFaces[eventI]
                << " sourceProc=" << acceptedSourceProcs[eventI]
                << " centre=" << acceptedCentres[eventI]
                << " seedRadius=" << acceptedRadii[eventI]
                << " stencilCells=" << stencilCells
                << " wallSuperheat=" << acceptedSuperheats[eventI]
                << " hazardRate=" << acceptedHazards[eventI]
                << " probability=" << acceptedProbabilities[eventI]
                << " couplingMode=" << sourceCouplingMode_
                << endl;
        }
    }

    return acceptedCentres.size();
}


void distributedSurfaceTopologyBirth::applyActiveEvents
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
)
{
    const fvMesh& mesh = phase1_.mesh();
    const scalar deltaT = mesh.time().deltaTValue();
    const scalar timeValue = mesh.time().value();
    const bool conservativeCoupling =
        sourceCouplingMode_ == "conservative";

    if (deltaT <= SMALL)
    {
        return;
    }

    forAll(seedCellLabels_, stencilI)
    {
        if (seedCellEndTimes_[stencilI] + SMALL < timeValue)
        {
            continue;
        }

        const label celli = seedCellLabels_[stencilI];
        const scalar targetVapor = seedCellTargetVapor_[stencilI];
        const scalar remainingDuration =
            max(seedCellEndTimes_[stencilI] - timeValue, deltaT);
        const scalar scheduleFraction =
            min(deltaT/remainingDuration, scalar(1));

        activeSeedMask_[celli] =
            max
            (
                activeSeedMask_[celli],
                targetVapor/max(targetSeedVaporFraction_, scalar(SMALL))
            );

        const scalar currentVapor =
            scalar(1) - phase1_[celli];
        const scalar vaporDeficit =
            targetVapor - currentVapor;

        if (vaporDeficit <= SMALL || phase1_[celli] <= SMALL)
        {
            continue;
        }

        scalar alphaStep =
            min
            (
                vaporDeficit*scheduleFraction,
                maxAlphaVaporPerStep_
            );
        alphaStep =
            min(alphaStep, max(phase1_[celli], scalar(0)));

        if (conservativeCoupling && thermalReserveFraction_ > 0)
        {
            const scalar latentHeat = satModel_.L()[celli];
            const scalar cellSuperheat =
                max
                (
                    liquidTemperature[celli] - TSat(celli),
                    scalar(0)
                );

            const scalar thermalAlphaLimit =
                latentHeat > SMALL
              ? thermalReserveFraction_
               *phase1_[celli]
               *liquidCp[celli]
               *cellSuperheat/latentHeat
              : scalar(0);

            alphaStep =
                min(alphaStep, max(thermalAlphaLimit, scalar(0)));
        }

        if (alphaStep <= SMALL)
        {
            continue;
        }

        alphaBirthSource_[celli] += alphaStep/deltaT;
    }

    integratedAlphaSource = 0;
    integratedMassSource = 0;
    integratedLatentSink = 0;
    integratedTopologyMassDefectRate = 0;
    integratedTopologyLatentDefectRate = 0;
    maxAlphaRate = 0;

    forAll(alphaBirthSource_, celli)
    {
        const scalar maxAlphaRateFromAvailableLiquid =
            max(phase1_[celli], scalar(0))/deltaT;

        if
        (
            alphaBirthSource_[celli]
          > maxAlphaRateFromAvailableLiquid
         && alphaBirthSource_[celli] > SMALL
        )
        {
            alphaBirthSource_[celli] =
                maxAlphaRateFromAvailableLiquid;
        }

        const scalar equivalentMassRate =
            liquidDensity[celli]*alphaBirthSource_[celli];
        const scalar equivalentLatentRate =
            equivalentMassRate*satModel_.L()[celli];

        if (conservativeCoupling)
        {
            massBirthSource_[celli] = equivalentMassRate;
            latentSink_[celli] = equivalentLatentRate;
        }
        else
        {
            equivalentMassDefectRate_[celli] =
                equivalentMassRate;
            equivalentLatentDefectRate_[celli] =
                equivalentLatentRate;
        }

        integratedAlphaSource +=
            alphaBirthSource_[celli]*mesh.V()[celli];
        integratedMassSource +=
            massBirthSource_[celli]*mesh.V()[celli];
        integratedLatentSink +=
            latentSink_[celli]*mesh.V()[celli];
        integratedTopologyMassDefectRate +=
            equivalentMassDefectRate_[celli]*mesh.V()[celli];
        integratedTopologyLatentDefectRate +=
            equivalentLatentDefectRate_[celli]*mesh.V()[celli];
        maxAlphaRate =
            max(maxAlphaRate, alphaBirthSource_[celli]);
    }
}


void distributedSurfaceTopologyBirth::updateSources()
{
    const fvMesh& mesh = phase1_.mesh();
    const label timeIndex = mesh.time().timeIndex();

    if (timeIndex == lastUpdateTimeIndex_)
    {
        return;
    }
    lastUpdateTimeIndex_ = timeIndex;

    superheat_ *= scalar(0);
    eligibleMask_ *= scalar(0);
    activeSeedMask_ *= scalar(0);
    alphaBirthSource_ *= scalar(0);
    massBirthSource_ *= scalar(0);
    latentSink_ *= scalar(0);
    equivalentMassDefectRate_ *= scalar(0);
    equivalentLatentDefectRate_ *= scalar(0);

    const volScalarField& liquidTemperature =
        phase1_.thermo().T();
    tmp<volScalarField> tLiquidDensity =
        phase1_.thermo().rho();
    const volScalarField& liquidDensity =
        tLiquidDensity();
    tmp<volScalarField> tLiquidCp =
        phase1_.thermo().Cp();
    const volScalarField& liquidCp =
        tLiquidCp();

    const label birthsThisStep =
        spawnEvents(liquidTemperature);

    scalar integratedAlphaSource = 0;
    scalar integratedMassSource = 0;
    scalar integratedLatentSink = 0;
    scalar integratedTopologyMassDefectRate = 0;
    scalar integratedTopologyLatentDefectRate = 0;
    scalar maxAlphaRate = 0;

    applyActiveEvents
    (
        liquidTemperature,
        liquidDensity,
        liquidCp,
        integratedAlphaSource,
        integratedMassSource,
        integratedLatentSink,
        integratedTopologyMassDefectRate,
        integratedTopologyLatentDefectRate,
        maxAlphaRate
    );

    // Convert local source integrals to global physical diagnostics. Each rank
    // receives the same reduced values, so cumulative defect state stays
    // identical across the decomposition.
    reduce(integratedAlphaSource, sumOp<scalar>());
    reduce(integratedMassSource, sumOp<scalar>());
    reduce(integratedLatentSink, sumOp<scalar>());
    reduce(integratedTopologyMassDefectRate, sumOp<scalar>());
    reduce(integratedTopologyLatentDefectRate, sumOp<scalar>());
    reduce(maxAlphaRate, maxOp<scalar>());

    const scalar deltaT = mesh.time().deltaTValue();
    const scalar topologyMassDefectThisStep =
        integratedTopologyMassDefectRate*deltaT;
    const scalar topologyLatentDefectThisStep =
        integratedTopologyLatentDefectRate*deltaT;

    cumulativeTopologyMassDefect_ +=
        topologyMassDefectThisStep;
    cumulativeTopologyLatentDefect_ +=
        topologyLatentDefectThisStep;

    superheat_.correctBoundaryConditions();
    eligibleMask_.correctBoundaryConditions();
    activeSeedMask_.correctBoundaryConditions();
    alphaBirthSource_.correctBoundaryConditions();
    massBirthSource_.correctBoundaryConditions();
    latentSink_.correctBoundaryConditions();
    equivalentMassDefectRate_.correctBoundaryConditions();
    equivalentLatentDefectRate_.correctBoundaryConditions();

    if (writeDiagnostics_)
    {
        label eligibleCells = 0;
        forAll(eligibleMask_, celli)
        {
            eligibleCells += eligibleMask_[celli] > 0.5;
        }
        reduce(eligibleCells, sumOp<label>());

        const auto maxWallSuperheat = gMax(superheat_);

        if (Pstream::master())
        {
            Info<< "DSTB"
                << " time=" << mesh.time().value()
                << " couplingMode=" << sourceCouplingMode_
                << " eligibleCells=" << eligibleCells
                << " births=" << birthsThisStep
                << " activeEvents="
                << activeEventCount(mesh.time().value())
                << " totalEvents=" << nextEventId_
                << " integratedAlphaSource="
                << integratedAlphaSource
                << " integratedMassSource="
                << integratedMassSource
                << " integratedLatentSink="
                << integratedLatentSink
                << " integratedTopologyMassDefectRate="
                << integratedTopologyMassDefectRate
                << " integratedTopologyLatentDefectRate="
                << integratedTopologyLatentDefectRate
                << " topologyMassDefectThisStep="
                << topologyMassDefectThisStep
                << " topologyLatentDefectThisStep="
                << topologyLatentDefectThisStep
                << " cumulativeTopologyMassDefect="
                << cumulativeTopologyMassDefect_
                << " cumulativeTopologyLatentDefect="
                << cumulativeTopologyLatentDefect_
                << " maxAlphaRate=" << maxAlphaRate
                << " maxWallSuperheat=" << maxWallSuperheat
                << endl;
        }
    }
}


void distributedSurfaceTopologyBirth::TSource1
(
    fvScalarMatrix& T1Eqn
)
{
    updateSources();
    T1Eqn.source() -=
        latentSink_.internalField()*phase1_.mesh().V();
}


void distributedSurfaceTopologyBirth::TSource2
(
    fvScalarMatrix& T2Eqn
)
{}


void distributedSurfaceTopologyBirth::energySource
(
    volScalarField& Q
)
{}


void distributedSurfaceTopologyBirth::energySource1
(
    volScalarField& q1
)
{}


void distributedSurfaceTopologyBirth::energySource2
(
    volScalarField& q2
)
{}


void distributedSurfaceTopologyBirth::massSource
(
    volScalarField& rhoSource
)
{
    updateSources();
    rhoSource += massBirthSource_;
}


void distributedSurfaceTopologyBirth::alphaSource
(
    volScalarField& rhoSource
)
{
    updateSources();
    rhoSource += alphaBirthSource_;
}

} // End namespace Foam

// ************************************************************************* //
