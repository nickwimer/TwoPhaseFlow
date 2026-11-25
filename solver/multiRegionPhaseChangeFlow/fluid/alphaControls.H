const dictionary& alphaControls = mesh.solverDict(alpha1.name());

label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));

// limitAlphaSource (default false): each alpha (sub)step takes no more liquid
// from a cell than the cell holds at its start, and adds no more than its
// vapor. See alphaSuSp.H.
const bool limitAlphaSource
(
    alphaControls.getOrDefault<bool>("limitAlphaSource", false)
);
