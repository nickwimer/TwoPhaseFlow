    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      - fvm::Sp(contErr, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    const dictionary& capillaryVerificationControlDict =
        mesh.time().controlDict();
    const Switch capillaryVerificationEnabled =
        capillaryVerificationControlDict.lookupOrDefault<Switch>
        (
            "capillaryVerificationEnabled",
            false
        );
    const scalar capillaryVerificationMomentumResidualScale =
        capillaryVerificationControlDict.lookupOrDefault<scalar>
        (
            "capillaryVerificationMomentumResidualScale",
            1.0
        );

    if (capillaryVerificationMomentumResidualScale < 0)
    {
        FatalErrorInFunction
            << "capillaryVerificationMomentumResidualScale must be >= 0: "
            << capillaryVerificationMomentumResidualScale
            << exit(FatalError);
    }

    surfaceScalarField capillaryFaceForceMomentumDiagnostic
    (
        "capillaryFaceForceMomentumDiagnostic",
        surfForces.surfaceTensionForce()
      + surfForces.accelerationForce()
    );
    surfaceScalarField pressureFaceGradientMomentumDiagnostic
    (
        "pressureFaceGradientMomentumDiagnostic",
        fvc::snGrad(p_rgh)
    );
    surfaceScalarField capillaryPressureResidualMomentumDiagnostic
    (
        "capillaryPressureResidualMomentumDiagnostic",
        capillaryFaceForceMomentumDiagnostic
      - pressureFaceGradientMomentumDiagnostic
    );

    scalar maxBoundaryFaceResidual = 0.0;
    const surfaceScalarField::Boundary& momentumResidualBoundary =
        capillaryPressureResidualMomentumDiagnostic.boundaryField();
    forAll(momentumResidualBoundary, patchi)
    {
        forAll(momentumResidualBoundary[patchi], facei)
        {
            maxBoundaryFaceResidual = max
            (
                maxBoundaryFaceResidual,
                mag(momentumResidualBoundary[patchi][facei])
            );
        }
    }
    reduce(maxBoundaryFaceResidual, maxOp<scalar>());

    tmp<volVectorField> tMomentumRhsDiagnostic = fvc::reconstruct
    (
        capillaryPressureResidualMomentumDiagnostic*mesh.magSf()
    );
    const volVectorField& momentumRhsDiagnostic = tMomentumRhsDiagnostic();

    tmp<volScalarField> tUEqnA = UEqn.A();
    const volScalarField& UEqnA = tUEqnA();
    tmp<volVectorField> tUEqnH = UEqn.H();
    const volVectorField& UEqnH = tUEqnH();

    if (capillaryVerificationEnabled)
    {
        Info<< "CAPMOMPRE time " << runTime.value()
            << " oCorr " << oCorr
            << " scale " << capillaryVerificationMomentumResidualScale
            << " maxAlphaDeltaOld "
            << gMax
               (
                   mag
                   (
                       alpha1.primitiveField()
                     - alpha1.oldTime().primitiveField()
                   )
               )
            << " maxContErr " << gMax(mag(contErr.field()))
            << " maxRhoPhi " << gMax(mag(rhoPhi.primitiveField()))
            << " maxFaceResidualInternal "
            << gMax
               (
                   mag
                   (
                       capillaryPressureResidualMomentumDiagnostic
                      .primitiveField()
                   )
               )
            << " maxFaceResidualBoundary " << maxBoundaryFaceResidual
            << " maxMomentumRhs "
            << gMax(mag(momentumRhsDiagnostic.primitiveField()))
            << " minA " << gMin(UEqnA.primitiveField())
            << " maxA " << gMax(UEqnA.primitiveField())
            << " maxH " << gMax(mag(UEqnH.primitiveField()))
            << " maxU " << gMax(mag(U.primitiveField()))
            << endl;
    }

    {
        const bool capillaryBalanceDiagnostics = pimple.lookupOrDefault
        (
            "capillaryBalanceDiagnostics",
            false
        );

        if (capillaryBalanceDiagnostics)
        {
            Info<< "CAPBAL preUEqn"
                << " time " << runTime.value()
                << " oCorr " << oCorr
                << " maxAlphaDeltaOld "
                << gMax
                   (
                       mag
                       (
                           alpha1.primitiveField()
                         - alpha1.oldTime().primitiveField()
                       )
                   )
                << " maxCapillary "
                << gMax(mag(capillaryFaceForceMomentumDiagnostic.primitiveField()))
                << " maxPressureGrad "
                << gMax
                   (
                       mag
                       (
                           pressureFaceGradientMomentumDiagnostic.primitiveField()
                       )
                   )
                << " maxFaceResidual "
                << gMax
                   (
                       mag
                       (
                           capillaryPressureResidualMomentumDiagnostic
                          .primitiveField()
                       )
                   )
                << " maxU " << gMax(mag(U.primitiveField()))
                << endl;
        }
    }

    if (momentumPredictor)
    {
        solve
        (
            UEqn
         == capillaryVerificationMomentumResidualScale*momentumRhsDiagnostic,
            mesh.solver(U.select(finalIter))
        );

        fvOptions.correct(U);
        K = 0.5*magSqr(U);

        if (capillaryVerificationEnabled)
        {
            Info<< "CAPMOMPOST time " << runTime.value()
                << " oCorr " << oCorr
                << " scale " << capillaryVerificationMomentumResidualScale
                << " maxU " << gMax(mag(U.primitiveField()))
                << endl;
        }
    }
