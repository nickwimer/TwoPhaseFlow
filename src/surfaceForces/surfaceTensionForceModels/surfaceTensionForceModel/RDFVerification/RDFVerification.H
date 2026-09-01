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

#ifndef RDFVerification_H
#define RDFVerification_H

#include "RDF.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class RDFVerification Declaration
\*---------------------------------------------------------------------------*/

// Verification-only wrapper around RDF.  It first computes the normal RDF
// curvature field and then optionally decomposes its interface curvature into
// mean bias and spatial-noise components.  Production cases must continue to
// select surfaceTensionForceModel RDF; the default RDF equations are not
// modified by this class.
class RDFVerification
:
    public RDF
{
    word mode_;
    scalar exactCurvature_;
    scalar alphaTol_;

    // Disallow copy construction/assignment
    RDFVerification(const RDFVerification&);
    void operator=(const RDFVerification&);

public:

    TypeName("RDFVerification");

    RDFVerification
    (
        const dictionary& dict,
        const volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U
    );

    virtual ~RDFVerification()
    {}

    virtual void correct();
};

} // End namespace Foam

#endif

// ************************************************************************* //
