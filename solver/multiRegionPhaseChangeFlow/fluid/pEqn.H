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

    if
    (
        capillaryVerificationVaporDensityTransportScale < 0
     || capillaryVerificationVaporPressureDdtScale < 0
     || capillaryVerificationVaporDensityUpdateScale < 0
    )
    {
        FatalErrorInFunction
            << "Case-202 vapor-coupling verification scales must be >= 0: "
            << capillaryVerificationVaporDensityTransportScale << ' '
            << capillaryVerificationVaporPressureDdtScale << ' '
            << capillaryVerificationVaporDensityUpdateScale
            << exit(FatalError);
    }

    if (corr == 0)
    {
        Info<< "CAPVFY vapor coupling scales: densityTransport="
            << capillaryVerificationVaporDensityTransportScale
            << " pressureDdt="
            << capillaryVerificationVaporPressureDdtScale
            << " densityUpdate="
            << capillaryVerificationVaporDensityUpdateScale
            << endl;
    }

    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi, UfFluid[i]))
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

            phi = phiHbyA + p_rghEqnIncomp.flux();

            U = HbyA
              + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);

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

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value()
        << "  average(p_rgh) " << p_rgh.weightedAverage(mesh.V()) << endl;

    if (mixture.thermo2().dpdt())
    {
       Info << "dpdt switch active" << endl;
       dpdt = fvc::ddt(p);
    }
}
