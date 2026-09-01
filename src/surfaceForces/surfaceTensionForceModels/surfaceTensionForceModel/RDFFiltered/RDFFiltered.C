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
#include "reconstructionSchemes.H"
#include "Switch.H"
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

    if (filterPasses_ > 0 && filterStrength_ > 0)
    {
        // RDF already extends K into a compact near-interface band. Use that
        // band as the filter support and explicitly include geometric cut
        // cells. The weighted neighbour average excludes zero-valued cells
        // outside the band, so a constant curvature field is preserved exactly
        // at band edges.
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

        // Use the filtered cell curvature in the same face interpolation and
        // CSF force path as RDF. No pressure, delta-function, or alpha equation
        // is changed by this model.
        Kf_ = fvc::interpolate(K_);
    }

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

    // Verification-only local oracle for an axis-aligned ellipsoid. The
    // reconstructed PLIC interface centre is used, so the comparison is made
    // at the location where RDF samples curvature rather than at the cell
    // centre. This diagnostic does not alter K, Kf, alpha, pressure, or force.
    const dictionary& controlDict = mesh.time().controlDict();
    const Switch ellipsoidDiagnostics =
        controlDict.lookupOrDefault<Switch>
        (
            "capillaryVerificationEllipsoidCurvatureDiagnostics",
            false
        );

    if (ellipsoidDiagnostics)
    {
        const vector origin =
            controlDict.lookupOrDefault<vector>
            (
                "capillaryVerificationEllipsoidOrigin",
                vector::zero
            );
        const vector semiAxis =
            controlDict.lookupOrDefault<vector>
            (
                "capillaryVerificationEllipsoidSemiAxis",
                vector(1, 1, 1)
            );

        if
        (
            semiAxis.x() <= VSMALL
         || semiAxis.y() <= VSMALL
         || semiAxis.z() <= VSMALL
        )
        {
            FatalErrorInFunction
                << "capillaryVerificationEllipsoidSemiAxis must be positive: "
                << semiAxis
                << exit(FatalError);
        }

        const reconstructionSchemes& surf =
            mesh.lookupObject<reconstructionSchemes>("reconstructionScheme");
        const volVectorField& interfaceCentre = surf.centre();
        const volVectorField& interfaceNormal = surf.normal();

        const scalar hx = 2.0/sqr(semiAxis.x());
        const scalar hy = 2.0/sqr(semiAxis.y());
        const scalar hz = 2.0/sqr(semiAxis.z());
        const scalar traceH = hx + hy + hz;

        label nOracle = 0;
        scalar sumOracle = 0.0;
        scalar sumOracleSq = 0.0;
        scalar sumK = 0.0;
        scalar sumKSq = 0.0;
        scalar sumKOracle = 0.0;
        scalar sumError = 0.0;
        scalar sumErrorSq = 0.0;
        scalar maxAbsError = 0.0;

        forAll(alpha, celli)
        {
            if
            (
                alpha[celli] > alphaTol_
             && alpha[celli] < 1.0 - alphaTol_
             && mag(interfaceNormal[celli]) > VSMALL
            )
            {
                const vector q = interfaceCentre[celli] - origin;
                const vector gradF
                (
                    hx*q.x(),
                    hy*q.y(),
                    hz*q.z()
                );
                const scalar gradMag = mag(gradF);

                if (gradMag > VSMALL)
                {
                    const scalar gradHGrad =
                        hx*sqr(gradF.x())
                      + hy*sqr(gradF.y())
                      + hz*sqr(gradF.z());
                    const scalar divNormal =
                        traceH/gradMag
                      - gradHGrad/(sqr(gradMag)*gradMag);
                    const scalar oracleK = -divNormal;
                    const scalar k = K_[celli];
                    const scalar error = k - oracleK;

                    ++nOracle;
                    sumOracle += oracleK;
                    sumOracleSq += sqr(oracleK);
                    sumK += k;
                    sumKSq += sqr(k);
                    sumKOracle += k*oracleK;
                    sumError += error;
                    sumErrorSq += sqr(error);
                    maxAbsError = max(maxAbsError, mag(error));
                }
            }
        }

        reduce(nOracle, sumOp<label>());
        reduce(sumOracle, sumOp<scalar>());
        reduce(sumOracleSq, sumOp<scalar>());
        reduce(sumK, sumOp<scalar>());
        reduce(sumKSq, sumOp<scalar>());
        reduce(sumKOracle, sumOp<scalar>());
        reduce(sumError, sumOp<scalar>());
        reduce(sumErrorSq, sumOp<scalar>());
        reduce(maxAbsError, maxOp<scalar>());

        if (!nOracle)
        {
            FatalErrorInFunction
                << "ellipsoid curvature oracle found no reconstructed "
                << "interface centres"
                << exit(FatalError);
        }

        const scalar invN = 1.0/scalar(nOracle);
        const scalar oracleMean = sumOracle*invN;
        const scalar kMean = sumK*invN;
        const scalar oracleVar =
            max(0.0, sumOracleSq*invN - sqr(oracleMean));
        const scalar kVar =
            max(0.0, sumKSq*invN - sqr(kMean));
        const scalar covariance =
            sumKOracle*invN - kMean*oracleMean;
        const scalar oracleStd = sqrt(oracleVar);
        const scalar kStd = sqrt(kVar);
        const scalar meanError = sumError*invN;
        const scalar rmsError = sqrt(sumErrorSq*invN);
        const scalar gain =
            oracleVar > VSMALL ? covariance/oracleVar : 0.0;
        const scalar correlation =
            (oracleVar > VSMALL && kVar > VSMALL)
          ? covariance/sqrt(oracleVar*kVar)
          : 0.0;

        Info<< "CAPKELLIPSOID time " << alpha1_.time().value()
            << " nInterface " << nOracle
            << " oracleMean " << oracleMean
            << " oracleStd " << oracleStd
            << " KMean " << kMean
            << " KStd " << kStd
            << " meanError " << meanError
            << " rmsError " << rmsError
            << " maxAbsError " << maxAbsError
            << " gain " << gain
            << " correlation " << correlation
            << endl;
    }
}


// ************************************************************************* //
