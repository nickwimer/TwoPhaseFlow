volScalarField::Internal Sp
(
    IOobject
    (
        "Sp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dgdt.dimensions(), Zero)
);

volScalarField::Internal Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dgdt.dimensions(), Zero)
);

const volScalarField& alphaSource = phaseChange.alphaSource();

// The phase-change alpha source is a rate set by phaseChange.correct() in the
// step before, and it is applied in full in each alpha subcycle. A macro model
// that converts a cell's last liquid sizes its rate from the liquid the cell
// held then; the advection of the first subcycle can carry some of that
// liquid out, and the second subcycle then takes more than is left, leaving
// alpha below zero (populationDetachedBubbleBirth at a bubble's completion).
// With limitAlphaSource each subcycle's source is capped by what the cell
// holds at its start, and the volume not converted in the subcycle is
// reported (ALPHA_SOURCE_LIMITED). Off, the source is applied in full
// wherever the cell has any liquid (or, for condensation, any vapor).
const scalar rDeltaTAlpha = 1.0/runTime.deltaTValue();
label nAlphaSourceLimited = 0;
scalar evaporationNotApplied = 0;   // [m3/s] of liquid the cells did not hold
scalar condensationNotApplied = 0;  // [m3/s] of vapor the cells did not hold

forAll(dgdt, celli)
{
    if (dgdt[celli] > 0.0)
    {
        Sp[celli] -= dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
        Su[celli] += dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
    }
    else if (dgdt[celli] < 0.0)
    {
        Sp[celli] += dgdt[celli]/max(alpha1[celli], 1e-4);
    }
    if (limitAlphaSource)
    {
        if (alphaSource[celli] > 0.0)
        {
            const scalar held = max(alpha1[celli], 0.0)*rDeltaTAlpha;
            if (alphaSource[celli] > held)
            {
                ++nAlphaSourceLimited;
                evaporationNotApplied +=
                    (alphaSource[celli] - held)*mesh.V()[celli];
            }
            Su[celli] -= min(alphaSource[celli], held);
        }
        else if (alphaSource[celli] < 0.0)
        {
            const scalar held = max(1.0 - alpha1[celli], 0.0)*rDeltaTAlpha;
            if (-alphaSource[celli] > held)
            {
                ++nAlphaSourceLimited;
                condensationNotApplied +=
                    (-alphaSource[celli] - held)*mesh.V()[celli];
            }
            Su[celli] += min(-alphaSource[celli], held);
        }
    }
    else if (alphaSource[celli] > 0.0 && alpha1[celli] > 0.0)
    {
        //Sp[celli] -= alphaSource[celli];
        Su[celli] -= alphaSource[celli];
    }
    else if (alphaSource[celli] < 0.0 && alpha1[celli] < 1.0)
    {
        Su[celli] -= alphaSource[celli];
    }

}

if (limitAlphaSource)
{
    reduce(nAlphaSourceLimited, sumOp<label>());
    if (nAlphaSourceLimited > 0)
    {
        reduce(evaporationNotApplied, sumOp<scalar>());
        reduce(condensationNotApplied, sumOp<scalar>());
        Info<< "ALPHA_SOURCE_LIMITED time=" << runTime.value()
            << " cells=" << nAlphaSourceLimited
            << " evaporationVolume="
            << evaporationNotApplied/rDeltaTAlpha
            << " condensationVolume="
            << condensationNotApplied/rDeltaTAlpha
            << endl;
    }
}

volScalarField::Internal divU
(
    mesh.moving()
  ? fvc::div(phi + mesh.phi())
  : fvc::div(phi)
);
