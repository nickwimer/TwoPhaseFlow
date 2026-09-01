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

#ifndef RDFFiltered_H
#define RDFFiltered_H

#include "RDF.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class RDFFiltered Declaration
\*---------------------------------------------------------------------------*/

// Experimental local curvature filter around the production RDF estimator.
// RDF::correct() is executed first.  A compact, convex neighbour filter is
// then applied only in the RDF near-interface curvature band.  The stencil is
// O(dx), preserves a constant curvature field exactly, and therefore vanishes
// as a physical smoothing length under mesh refinement.  Production cases
// continue to select surfaceTensionForceModel RDF until this model is fully
// verified on static and dynamic interface benchmarks.
class RDFFiltered
:
    public RDF
{
    scalar filterStrength_;
    label filterPasses_;
    scalar alphaTol_;

    // Disallow copy construction/assignment
    RDFFiltered(const RDFFiltered&);
    void operator=(const RDFFiltered&);

public:

    TypeName("RDFFiltered");

    RDFFiltered
    (
        const dictionary& dict,
        const volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U
    );

    virtual ~RDFFiltered()
    {}

    virtual void correct();
};

} // End namespace Foam

#endif

// ************************************************************************* //
