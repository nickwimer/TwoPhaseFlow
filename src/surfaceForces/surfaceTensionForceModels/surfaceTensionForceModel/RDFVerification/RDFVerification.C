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

#include "RDFVerification.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"

namespace Foam
{
    defineTypeNameAndDebug(RDFVerification, 0);
    addToRunTimeSelectionTable
    (
        surfaceTensionForceModel,
        RDFVerification,
        components
    );
}


Foam::RDFVerification::RDFVerification
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    RDF(dict, alpha1, phi, U),
    mode_(dict.lookupOrDefault<word>("rdfVerificationMode", "raw")),
    exactCurvature_
    (
        dict.lookupOrDefault<scalar>
        (
            "rdfVerificationExactCurvature",
            -4000.0
        )
    ),
    alphaTol_
    (
        dict.lookupOrDefault<scalar>
        (
            "rdfVerificationAlphaTol",
            1e-8
        )
    )
{
    if
    (
        mode_ != "raw"
     && mode_ != "meanOnly"
     && mode_ != "debiased"
     && mode_ != "halfNoise"
    )
    {
        FatalErrorInFunction
            << "Unknown rdfVerificationMode '" << mode_ << "'. Valid modes: "
            << "raw meanOnly debiased halfNoise"
            << exit(FatalError);
    }

    if (alphaTol_ <= 0 || alphaTol_ >= 0.5)
    {
        FatalErrorInFunction
            << "rdfVerificationAlphaTol must be in (0, 0.5): "
            << alphaTol_
            << exit(FatalError);
    }
}


void Foam::RDFVerification::correct()
{
    // First execute the unmodified production RDF algorithm.
    RDF::correct();

    const scalarField& alpha = alpha1_.primitiveField();
    const scalarField& rawCurvature = K_.primitiveField();

    label nInterface = 0;
    scalar sumCurvature = 0.0;

    forAll(alpha, celli)
    {
        if (alpha[celli] > alphaTol_ && alpha[celli] < 1.0 - alphaTol_)
        {
            ++nInterface;
            sumCurvature += rawCurvature[celli];
        }
    }

    reduce(nInterface, sumOp<label>());
    reduce(sumCurvature, sumOp<scalar>());

    if (!nInterface)
    {
        FatalErrorInFunction
            << "RDFVerification found no interface cells using alphaTol="
            << alphaTol_
            << exit(FatalError);
    }

    const scalar rawMean = sumCurvature/scalar(nInterface);

    if (mode_ == "raw")
    {
        Info<< "CAPKXFORM time " << alpha1_.time().value()
            << " mode " << mode_
            << " rawMean " << rawMean
            << " targetMean " << rawMean
            << " noiseScale 1"
            << " nInterface " << nInterface
            << endl;
        return;
    }

    scalar targetMean = exactCurvature_;
    scalar noiseScale = 1.0;

    if (mode_ == "meanOnly")
    {
        targetMean = rawMean;
        noiseScale = 0.0;
    }
    else if (mode_ == "debiased")
    {
        targetMean = exactCurvature_;
        noiseScale = 1.0;
    }
    else if (mode_ == "halfNoise")
    {
        targetMean = exactCurvature_;
        noiseScale = 0.5;
    }

    // Keep the cell curvature diagnostics consistent with the transformed
    // field.  Only cells where RDF supplied curvature are altered; cells
    // outside the RDF near-interface band remain zero.
    scalarField& curvature = K_.primitiveFieldRef();
    forAll(curvature, celli)
    {
        if (mag(curvature[celli]) > VSMALL)
        {
            curvature[celli] =
                targetMean + noiseScale*(curvature[celli] - rawMean);
        }
    }

    // Transform the actual face curvature used by sigma*Kf*snGrad(alpha).
    // Restrict the operation to force-active faces so values away from the
    // interface remain untouched.
    const surfaceScalarField& delta = deltaFunctionModel_->deltaFunction();
    scalarField& curvatureFaces = Kf_.primitiveFieldRef();
    const scalarField& deltaFaces = delta.primitiveField();

    forAll(curvatureFaces, facei)
    {
        if (mag(deltaFaces[facei]) > VSMALL)
        {
            curvatureFaces[facei] =
                targetMean + noiseScale*(curvatureFaces[facei] - rawMean);
        }
    }

    surfaceScalarField::Boundary& curvatureBoundary = Kf_.boundaryFieldRef();
    const surfaceScalarField::Boundary& deltaBoundary = delta.boundaryField();

    forAll(curvatureBoundary, patchi)
    {
        forAll(curvatureBoundary[patchi], facei)
        {
            if (mag(deltaBoundary[patchi][facei]) > VSMALL)
            {
                curvatureBoundary[patchi][facei] =
                    targetMean
                  + noiseScale
                   *(
                        curvatureBoundary[patchi][facei]
                      - rawMean
                    );
            }
        }
    }

    Info<< "CAPKXFORM time " << alpha1_.time().value()
        << " mode " << mode_
        << " rawMean " << rawMean
        << " targetMean " << targetMean
        << " noiseScale " << noiseScale
        << " nInterface " << nInterface
        << endl;
}


// ************************************************************************* //
