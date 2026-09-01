/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
    unofficial extension to OpenFOAM.
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
\*---------------------------------------------------------------------------*/

#include "RDFFiltered.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"
#include "fvc.H"

namespace Foam
{
    defineTypeNameAndDebug(RDFFiltered, 0);
    addToRunTimeSelectionTable
    (
        surfaceTensionForceModel,
        RDFFiltered,
        components
    );
}


Foam::RDFFiltered::RDFFiltered
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    RDF(dict, alpha1, phi, U),
    filterStrength_
    (
        dict.lookupOrDefault<scalar>("curvatureFilterStrength", 0.5)
    ),
    filterPasses_
    (
        dict.lookupOrDefault<label>("curvatureFilterPasses", 1)
    ),
    alphaTol_
    (
        dict.lookupOrDefault<scalar>("curvatureFilterAlphaTol", 1e-8)
    )
{
    if (filterStrength_ < 0 || filterStrength_ > 1)
    {
        FatalErrorInFunction
            << "curvatureFilterStrength must be in [0, 1]: "
            << filterStrength_
            << exit(FatalError);
    }

    if (filterPasses_ < 0)
    {
        FatalErrorInFunction
            << "curvatureFilterPasses must be >= 0: "
            << filterPasses_
            << exit(FatalError);
    }

    if (alphaTol_ <= 0 || alphaTol_ >= 0.5)
    {
        FatalErrorInFunction
            << "curvatureFilterAlphaTol must be in (0, 0.5): "
            << alphaTol_
            << exit(FatalError);
    }
}


void Foam::RDFFiltered::correct()
{
    // Always start from the unmodified production RDF curvature.
    RDF::correct();

    const fvMesh& mesh = alpha1_.mesh();
    const scalarField& alpha = alpha1_.primitiveField();

    label nInterface = 0;
    scalar rawSum = 0.0;
    scalar rawSumSq = 0.0;

    forAll(alpha, celli)
    {
        if (alpha[celli] > alphaTol_ && alpha[celli] < 1.0 - alphaTol_)
        {
            const scalar k = K_[celli];
            ++nInterface;
            rawSum += k;
            rawSumSq += k*k;
        }
    }

    reduce(nInterface, sumOp<label>());
    reduce(rawSum, sumOp<scalar>());
    reduce(rawSumSq, sumOp<scalar>());

    if (!nInterface)
    {
        FatalErrorInFunction
            << "RDFFiltered found no interface cells using alphaTol="
            << alphaTol_
            << exit(FatalError);
    }

    const scalar rawMean = rawSum/scalar(nInterface);
    const scalar rawVariance =
        max(0.0, rawSumSq/scalar(nInterface) - sqr(rawMean));
    const scalar rawStd = sqrt(rawVariance);

    if (filterPasses_ == 0 || filterStrength_ == 0)
    {
        Info<< "CAPKFILTER time " << alpha1_.time().value()
            << " strength " << filterStrength_
            << " passes " << filterPasses_
            << " nInterface " << nInterface
            << " rawMean " << rawMean
            << " rawStd " << rawStd
            << " filteredMean " << rawMean
            << " filteredStd " << rawStd
            << endl;
        return;
    }

    // RDF already extends K into a compact near-interface band.  Use that band
    // as the filter support and explicitly include geometric cut cells.  The
    // weighted neighbour average below excludes zero-valued cells outside the
    // band, so a constant curvature field is preserved exactly at band edges.
    volScalarField activeMask
    (
        pos0
        (
            mag(K_)
          - dimensionedScalar("kActive", K_.dimensions(), VSMALL)
        )
    );

    scalarField& active = activeMask.primitiveFieldRef();
    forAll(alpha, celli)
    {
        if (alpha[celli] > alphaTol_ && alpha[celli] < 1.0 - alphaTol_)
        {
            active[celli] = 1.0;
        }
    }
    activeMask.correctBoundaryConditions();

    for (label pass = 0; pass < filterPasses_; ++pass)
    {
        volScalarField weightedCurvature
        (
            IOobject
            (
                "rdfFilterWeightedCurvature",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            activeMask*K_
        );

        tmp<surfaceScalarField> tMaskf = fvc::interpolate(activeMask);
        tmp<surfaceScalarField> tWeightedKf =
            fvc::interpolate(weightedCurvature);
        tmp<volScalarField> tDenominator = fvc::average(tMaskf());
        tmp<volScalarField> tNumerator = fvc::average(tWeightedKf());

        const scalarField oldCurvature(K_.primitiveField());
        scalarField& curvature = K_.primitiveFieldRef();
        const scalarField& denominator = tDenominator().primitiveField();
        const scalarField& numerator = tNumerator().primitiveField();

        forAll(curvature, celli)
        {
            if (active[celli] > 0.5 && denominator[celli] > SMALL)
            {
                const scalar neighbourMean =
                    numerator[celli]/denominator[celli];
                curvature[celli] =
                    (1.0 - filterStrength_)*oldCurvature[celli]
                  + filterStrength_*neighbourMean;
            }
        }

        K_.correctBoundaryConditions();
    }

    // Use the filtered cell curvature in the same face interpolation and CSF
    // force path as RDF.  No pressure, delta-function, or alpha equation is
    // changed by this model.
    Kf_ = fvc::interpolate(K_);

    scalar filteredSum = 0.0;
    scalar filteredSumSq = 0.0;

    forAll(alpha, celli)
    {
        if (alpha[celli] > alphaTol_ && alpha[celli] < 1.0 - alphaTol_)
        {
            const scalar k = K_[celli];
            filteredSum += k;
            filteredSumSq += k*k;
        }
    }

    reduce(filteredSum, sumOp<scalar>());
    reduce(filteredSumSq, sumOp<scalar>());

    const scalar filteredMean = filteredSum/scalar(nInterface);
    const scalar filteredVariance =
        max(0.0, filteredSumSq/scalar(nInterface) - sqr(filteredMean));
    const scalar filteredStd = sqrt(filteredVariance);

    Info<< "CAPKFILTER time " << alpha1_.time().value()
        << " strength " << filterStrength_
        << " passes " << filterPasses_
        << " nInterface " << nInterface
        << " rawMean " << rawMean
        << " rawStd " << rawStd
        << " filteredMean " << filteredMean
        << " filteredStd " << filteredStd
        << endl;
}


// ************************************************************************* //
