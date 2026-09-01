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

    if (capillaryVerificationEnabled)
    {
        Info<< "CAPMOMPRE time " << runTime.value()
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
            << " maxContErr " << gMax(mag(contErr))
            << " maxRhoPhi " << gMax(mag(rhoPhi.primitiveField()))
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
         ==
            fvc::reconstruct
            (
                (
                    //  mixture.surfaceTensionForce()
                    surfForces.surfaceTensionForce()
                    //- ghf*fvc::snGrad(rho)
                    + surfForces.accelerationForce()
                    - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            ),
            mesh.solver(U.select(finalIter))
        );

        fvOptions.correct(U);
        K = 0.5*magSqr(U);

        if (capillaryVerificationEnabled)
        {
            Info<< "CAPMOMPOST time " << runTime.value()
                << " oCorr " << oCorr
                << " maxU " << gMax(mag(U.primitiveField()))
                << endl;
        }
    }
