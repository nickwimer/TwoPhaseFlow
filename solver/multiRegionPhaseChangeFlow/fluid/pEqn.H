{

    bool closedVolume = p_rgh.needReference();
    dimensionedScalar compressibility = fvc::domainIntegrate(psi);
    bool compressible = (compressibility.value() > SMALL);
    const bool capillaryBalanceDiagnostics = pimple.lookupOrDefault
    (
        "capillaryBalanceDiagnostics",
        false
    );

    // Verification-only controls used by Case 202 to isolate the inherited
    // compressible-VOF pressure coupling. Read them from the top-level
    // controlDict so their scope is unambiguous in this multi-region solver.
    // Defaults preserve the production equations exactly. These are not
    // physical model switches and must remain at 1 outside verification.
    const dictionary& verificationControlDict = mesh.time().controlDict();
    const Switch capillaryVerificationEnabled =
        verificationControlDict.lookupOrDefault<Switch>
        (
            "capillaryVerificationEnabled",
            false
        );
    const scalar capillaryVerificationVaporDensityTransportScale =
        verificationControlDict.lookupOrDefault<scalar>
        (
            "capillaryVerificationVaporDensityTransportScale",
            1.0
        );
    const scalar capillaryVerificationVaporPressureDdtScale =
        verificationControlDict.lookupOrDefault<scalar>
        (
            "capillaryVerificationVaporPressureDdtScale",
            1.0
        );
    const scalar capillaryVerificationVaporDensityUpdateScale =
        verificationControlDict.lookupOrDefault<scalar>
        (
            "capillaryVerificationVaporDensityUpdateScale",
            1.0
        );
    const scalar capillaryVerificationDgdtScale =
        verificationControlDict.lookupOrDefault<scalar>
        (
            "capillaryVerificationDgdtScale",
            1.0
        );
    const scalar capillaryVerificationDdtCorrScale =
        verificationControlDict.lookupOrDefault<scalar>
        (
            "capillaryVerificationDdtCorrScale",
            1.0
        );

    if
    (
        capillaryVerificationVaporDensityTransportScale < 0
     || capillaryVerificationVaporPressureDdtScale < 0
     || capillaryVerificationVaporDensityUpdateScale < 0
     || capillaryVerificationDgdtScale < 0
     || capillaryVerificationDdtCorrScale < 0
    )
    {
        FatalErrorInFunction
            << "Case-202 verification scales must be >= 0: "
            << capillaryVerificationVaporDensityTransportScale << ' '
            << capillaryVerificationVaporPressureDdtScale << ' '
            << capillaryVerificationVaporDensityUpdateScale << ' '
            << capillaryVerificationDgdtScale << ' '
            << capillaryVerificationDdtCorrScale
            << exit(FatalError);
    }

    if (capillaryVerificationEnabled && corr == 0)
    {
        Info<< "CAPVFY vapor coupling scales: densityTransport="
            << capillaryVerificationVaporDensityTransportScale
            << " pressureDdt="
            << capillaryVerificationVaporPressureDdtScale
            << " densityUpdate="
            << capillaryVerificationVaporDensityUpdateScale
            << " dgdt=" << capillaryVerificationDgdtScale
            << " ddtCorr=" << capillaryVerificationDdtCorrScale
            << endl;
    }

    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField rawDdtCorrFlux
    (
        "rawDdtCorrFlux",
        MRF.zeroFilter
        (
            fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi, UfFluid[i])
        )
    );

    if (capillaryVerificationEnabled)
    {
        Info<< "CAPDDTCORR time " << runTime.value()
            << " oCorr " << oCorr
            << " corr " << corr
            << " scale " << capillaryVerificationDdtCorrScale
            << " rawMaxAbs "
            << gMax(mag(rawDdtCorrFlux.primitiveField()))
            << " scaledMaxAbs "
            << capillaryVerificationDdtCorrScale
               *gMax(mag(rawDdtCorrFlux.primitiveField()))
            << endl;
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + capillaryVerificationDdtCorrScale*rawDdtCorrFlux
    );
    MRF.makeRelative(phiHbyA);

    surfaceScalarField phig
    (
        (
            surfForces.surfaceTensionForce()
            + surfForces.accelerationForce()
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phiHbyA, U);

    tmp<fvScalarMatrix> p_rghEqnComp1;
    tmp<fvScalarMatrix> p_rghEqnComp2;

    {
        #include "rhofs.H"

        p_rghEqnComp1 =
            pos(alpha1)
           *(
                (
                    fvc::ddt(alpha1, rho1) + fvc::div(alphaPhi1*rho1f)
                  - (fvOptions(alpha1, mixture.thermo1().rho())&rho1)
                )/rho1
              - fvc::ddt(alpha1) - fvc::div(alphaPhi1)
              + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh))
            );

        p_rghEqnComp2 =
            pos(alpha2)
           *(
                capillaryVerificationVaporDensityTransportScale
               *(
                    (
                        fvc::ddt(alpha2, rho2) + fvc::div(alphaPhi2*rho2f)
                      - (fvOptions(alpha2, mixture.thermo2().rho())&rho2)
                    )/rho2
                  - fvc::ddt(alpha2) - fvc::div(alphaPhi2)
                )
              + capillaryVerificationVaporPressureDdtScale
               *(alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh))
            );
    }

    if (mesh.moving())
    {
        p_rghEqnComp1.ref() += fvc::div(mesh.phi())*alpha1;
        p_rghEqnComp2.ref() += fvc::div(mesh.phi())*alpha2;
    }

    p_rghEqnComp1.ref() *= pos(alpha1);
    p_rghEqnComp2.ref() *= pos(alpha2);

    // Cache p_rgh prior to solve for density update and diagnostics
    volScalarField p_rgh_0(p_rgh);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
            - fvm::laplacian(rAUf, p_rgh)
            - phaseChange.massSource()/rho
        );

        solve
        (
            p_rghEqnComp1() + p_rghEqnComp2() + p_rghEqnIncomp,
            mesh.solver
            (
                p_rgh.select
                (
                    oCorr == nOuterCorr-1
                    && corr == nCorr-1
                    && nonOrth == nNonOrthCorr
                )
            )
        );

        if (nonOrth == nNonOrthCorr)
        {
            p = max(p_rgh + (alpha1*rho1 + alpha2*rho2)*gh, pMin);
            p_rgh = p - (alpha1*rho1 + alpha2*rho2)*gh;

            dgdt =
            (
                alpha1*(p_rghEqnComp2 & p_rgh)
              - alpha2*(p_rghEqnComp1 & p_rgh)
            );

            const scalar rawDgdtMin = gMin(dgdt.primitiveField());
            const scalar rawDgdtMax = gMax(dgdt.primitiveField());
            const scalar rawDgdtMaxAbs = gMax(mag(dgdt.primitiveField()));

            dgdt *= capillaryVerificationDgdtScale;

            if (capillaryVerificationEnabled)
            {
                Info<< "CAPDGDT time " << runTime.value()
                    << " oCorr " << oCorr
                    << " corr " << corr
                    << " scale " << capillaryVerificationDgdtScale
                    << " rawMin " << rawDgdtMin
                    << " rawMax " << rawDgdtMax
                    << " rawMaxAbs " << rawDgdtMaxAbs
                    << " scaledMaxAbs "
                    << gMax(mag(dgdt.primitiveField()))
                    << endl;
            }

            phi = phiHbyA + p_rghEqnIncomp.flux();

            U = HbyA
              + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);

            if (capillaryVerificationEnabled)
            {
                surfaceScalarField phiUMismatch
                (
                    "phiUMismatch",
                    phi - fvc::flux(U)
                );
                Info<< "CAPPHIU time " << runTime.value()
                    << " oCorr " << oCorr
                    << " corr " << corr
                    << " maxAbs "
                    << gMax(mag(phiUMismatch.primitiveField()))
                    << endl;
            }

            if (capillaryBalanceDiagnostics)
            {
                surfaceScalarField capillaryFaceForce
                (
                    "capillaryFaceForcePostPDiagnostic",
                    surfForces.surfaceTensionForce()
                  + surfForces.accelerationForce()
                );
                surfaceScalarField pressureFaceGradient
                (
                    "pressureFaceGradientPostPDiagnostic",
                    fvc::snGrad(p_rgh)
                );
                surfaceScalarField capillaryPressureResidual
                (
                    "capillaryPressureResidualPostPDiagnostic",
                    capillaryFaceForce - pressureFaceGradient
                );

                Info<< "CAPBAL postPEqn"
                    << " time " << runTime.value()
                    << " oCorr " << oCorr
                    << " corr " << corr
                    << " maxPDelta "
                    << gMax
                       (
                           mag
                           (
                               p_rgh.primitiveField()
                             - p_rgh_0.primitiveField()
                           )
                       )
                    << " maxCapillary "
                    << gMax(mag(capillaryFaceForce.primitiveField()))
                    << " maxPressureGrad "
                    << gMax(mag(pressureFaceGradient.primitiveField()))
                    << " maxFaceResidual "
                    << gMax(mag(capillaryPressureResidual.primitiveField()))
                    << " maxU " << gMax(mag(U.primitiveField()))
                    << endl;
            }
        }
    }

    // Correct Uf if the mesh is moving
    fvc::correctUf(UfFluid[i], U, fvc::absolute(phi, U));

    if (closedVolume && compressible)
    {
        dimensionedScalar incompMass = initialMass*0.0;
        if(mixture.thermo1().incompressible())
        {
            incompMass = fvc::domainIntegrate(alpha1*rho1);
        }
        p +=
        (
            initialMass
            - incompMass
            - fvc::domainIntegrate(psi*p)
        )/compressibility;
        mixture.correctThermo(TL,TV);
        mixture.correct();
        rho = alpha1*rho1 + alpha2*rho2;
        p_rgh = p - rho*gh;
        p_rgh.correctBoundaryConditions();
    }
    else
    {
        // Update densities from change in p_rgh. The vapor scale is a
        // verification-only Case-202 control; production behavior is 1.0.
        mixture.thermo1().correctRho(psi1*(p_rgh - p_rgh_0));
        mixture.thermo2().correctRho
        (
            capillaryVerificationVaporDensityUpdateScale
           *psi2*(p_rgh - p_rgh_0)
        );
    }

    rho = alpha1*rho1 + alpha2*rho2;

    K = 0.5*magSqr(U);

    // Correct p_rgh for consistency with p and the updated densities
    p = max(p_rgh + rho*gh, pMin);
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();

    if (capillaryVerificationEnabled)
    {
        surfaceScalarField capillaryFaceForcePressurePatchDiagnostic
        (
            "capillaryFaceForcePressurePatchDiagnostic",
            surfForces.surfaceTensionForce()
          + surfForces.accelerationForce()
        );
        surfaceScalarField pressureFaceGradientPressurePatchDiagnostic
        (
            "pressureFaceGradientPressurePatchDiagnostic",
            fvc::snGrad(p_rgh)
        );
        surfaceScalarField capillaryPressureResidualPressurePatchDiagnostic
        (
            "capillaryPressureResidualPressurePatchDiagnostic",
            capillaryFaceForcePressurePatchDiagnostic
          - pressureFaceGradientPressurePatchDiagnostic
        );

        const surfaceScalarField::Boundary& capillaryBoundary =
            capillaryFaceForcePressurePatchDiagnostic.boundaryField();
        const surfaceScalarField::Boundary& pressureGradientBoundary =
            pressureFaceGradientPressurePatchDiagnostic.boundaryField();
        const surfaceScalarField::Boundary& residualBoundary =
            capillaryPressureResidualPressurePatchDiagnostic.boundaryField();
        const volScalarField::Boundary& pRghBoundary = p_rgh.boundaryField();

        forAll(pRghBoundary, patchi)
        {
            label nPatchFaces = pRghBoundary[patchi].size();
            reduce(nPatchFaces, sumOp<label>());
            if (!nPatchFaces)
            {
                continue;
            }

            scalar pBoundaryMin = GREAT;
            scalar pBoundaryMax = -GREAT;
            scalar pInternalMin = GREAT;
            scalar pInternalMax = -GREAT;
            scalar maxDeltaP = 0.0;
            scalar maxCapillary = 0.0;
            scalar maxPressureGradient = 0.0;
            scalar maxResidual = 0.0;
            const labelUList& faceCells = mesh.boundary()[patchi].faceCells();

            forAll(pRghBoundary[patchi], facei)
            {
                const scalar pBoundaryValue = pRghBoundary[patchi][facei];
                const scalar pInternalValue = p_rgh[faceCells[facei]];
                pBoundaryMin = min(pBoundaryMin, pBoundaryValue);
                pBoundaryMax = max(pBoundaryMax, pBoundaryValue);
                pInternalMin = min(pInternalMin, pInternalValue);
                pInternalMax = max(pInternalMax, pInternalValue);
                maxDeltaP = max
                (
                    maxDeltaP,
                    mag(pBoundaryValue - pInternalValue)
                );
                maxCapillary = max
                (
                    maxCapillary,
                    mag(capillaryBoundary[patchi][facei])
                );
                maxPressureGradient = max
                (
                    maxPressureGradient,
                    mag(pressureGradientBoundary[patchi][facei])
                );
                maxResidual = max
                (
                    maxResidual,
                    mag(residualBoundary[patchi][facei])
                );
            }

            reduce(pBoundaryMin, minOp<scalar>());
            reduce(pBoundaryMax, maxOp<scalar>());
            reduce(pInternalMin, minOp<scalar>());
            reduce(pInternalMax, maxOp<scalar>());
            reduce(maxDeltaP, maxOp<scalar>());
            reduce(maxCapillary, maxOp<scalar>());
            reduce(maxPressureGradient, maxOp<scalar>());
            reduce(maxResidual, maxOp<scalar>());

            Info<< "CAPPATCHP time " << runTime.value()
                << " oCorr " << oCorr
                << " corr " << corr
                << " patch " << mesh.boundary()[patchi].name()
                << " pType " << pRghBoundary[patchi].type()
                << " pBoundaryMin " << pBoundaryMin
                << " pBoundaryMax " << pBoundaryMax
                << " pInternalMin " << pInternalMin
                << " pInternalMax " << pInternalMax
                << " maxDeltaP " << maxDeltaP
                << " maxCapillary " << maxCapillary
                << " maxPressureGrad " << maxPressureGradient
                << " maxResidual " << maxResidual
                << endl;
        }
    }

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value()
        << "  average(p_rgh) " << p_rgh.weightedAverage(mesh.V()) << endl;

    if (mixture.thermo2().dpdt())
    {
       Info << "dpdt switch active" << endl;
       dpdt = fvc::ddt(p);
    }
}
