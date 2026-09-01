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
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::surfaceForces

Description
    computes the forces accting at the interace e.g. surface tension and
    gravity

SourceFiles
    surfaceForces.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceForces_H
#define surfaceForces_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"
#include "Switch.H"
#include "PstreamReduceOps.H"

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "surfaceTensionModel.H"
#include "surfaceTensionForceModel.H"
#include "deltaFunctionModel.H"
#include "accelerationForceModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class surfaceForces Declaration
\*---------------------------------------------------------------------------*/

class surfaceForces
{
private:

    //-
        dictionary surfaceForcesCoeffs_;

        //- reference to volume fraction field
        const volScalarField& alpha1_;

        const fvMesh& mesh_;


        autoPtr<surfaceTensionForceModel> surfTenForceModel_;

        autoPtr<accelerationForceModel> accModel_;


    // Private Member Functions

        //- Disallow copy construct
        surfaceForces(const surfaceForces&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceForces&);

public:


    // Constructors

        //- Construct from components
        surfaceForces
        (
            const volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfaceForces()
    {}


    // Member Functions
        virtual scalar capillaryDt
        (
            const volScalarField& rho1,
            const volScalarField& rho2
        );

        virtual scalar capillaryDt
        (
            const dimensionedScalar rho1,
            const dimensionedScalar rho2
        );

        dictionary& modelDict();

        const dictionary& modelDict() const;

        void correct()
        {
            surfTenForceModel_->correct();
            accModel_->correct();

            // Verification-only curvature-oracle diagnostics for canonical
            // interfaces such as Case 202. Defaults are disabled and do not
            // alter the selected surface-tension model or its force.
            const dictionary& controlDict = mesh_.time().controlDict();
            const Switch curvatureDiagnostics =
                controlDict.lookupOrDefault<Switch>
                (
                    "capillaryVerificationCurvatureDiagnostics",
                    false
                );

            if (curvatureDiagnostics)
            {
                const scalar exactCurvature =
                    controlDict.lookupOrDefault<scalar>
                    (
                        "capillaryVerificationExactCurvature",
                        0.0
                    );
                const scalar alphaTol =
                    controlDict.lookupOrDefault<scalar>
                    (
                        "capillaryVerificationInterfaceAlphaTol",
                        1e-8
                    );

                if (alphaTol <= 0 || alphaTol >= 0.5)
                {
                    FatalErrorInFunction
                        << "capillaryVerificationInterfaceAlphaTol must be "
                        << "in (0, 0.5): " << alphaTol
                        << exit(FatalError);
                }

                const volScalarField& curvature = surfTenForceModel_->K();
                label nInterface = 0;
                scalar minCurvature = GREAT;
                scalar maxCurvature = -GREAT;
                scalar sumCurvature = 0.0;
                scalar sumAbsError = 0.0;
                scalar sumSqError = 0.0;
                scalar maxAbsError = 0.0;

                const scalarField& alphaInternal = alpha1_.primitiveField();
                const scalarField& curvatureInternal = curvature.primitiveField();

                forAll(alphaInternal, celli)
                {
                    const scalar alpha = alphaInternal[celli];
                    if (alpha > alphaTol && alpha < 1.0 - alphaTol)
                    {
                        const scalar k = curvatureInternal[celli];
                        const scalar error = k - exactCurvature;
                        ++nInterface;
                        minCurvature = min(minCurvature, k);
                        maxCurvature = max(maxCurvature, k);
                        sumCurvature += k;
                        sumAbsError += mag(error);
                        sumSqError += error*error;
                        maxAbsError = max(maxAbsError, mag(error));
                    }
                }

                reduce(nInterface, sumOp<label>());
                reduce(minCurvature, minOp<scalar>());
                reduce(maxCurvature, maxOp<scalar>());
                reduce(sumCurvature, sumOp<scalar>());
                reduce(sumAbsError, sumOp<scalar>());
                reduce(sumSqError, sumOp<scalar>());
                reduce(maxAbsError, maxOp<scalar>());

                if (nInterface > 0)
                {
                    const scalar meanCurvature =
                        sumCurvature/scalar(nInterface);
                    const scalar meanAbsError =
                        sumAbsError/scalar(nInterface);
                    const scalar rmsError =
                        sqrt(sumSqError/scalar(nInterface));

                    Info<< "CAPK time " << mesh_.time().value()
                        << " exactK " << exactCurvature
                        << " alphaTol " << alphaTol
                        << " nInterface " << nInterface
                        << " minK " << minCurvature
                        << " maxK " << maxCurvature
                        << " meanK " << meanCurvature
                        << " meanAbsError " << meanAbsError
                        << " rmsError " << rmsError
                        << " maxAbsError " << maxAbsError
                        << endl;
                }
                else
                {
                    Info<< "CAPK time " << mesh_.time().value()
                        << " exactK " << exactCurvature
                        << " alphaTol " << alphaTol
                        << " nInterface 0"
                        << endl;
                }
            }
        }

        const volScalarField& acc() const
        {
            return accModel_->acc();
        }

        const surfaceScalarField& accf() const
        {
            return accModel_->accf();
        }

        tmp<volScalarField> sigma() const
        {
            return surfTenForceModel_->sigma();
        }

        const surfaceScalarField& nHatf() const
        {
            return surfTenForceModel_->nHatf();
        }

        tmp<surfaceScalarField> surfaceTensionForce();

        tmp<surfaceScalarField> accelerationForce();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
