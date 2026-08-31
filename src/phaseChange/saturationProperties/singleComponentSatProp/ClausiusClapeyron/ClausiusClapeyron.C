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

\*---------------------------------------------------------------------------*/

#include "ClausiusClapeyron.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

namespace Foam
{
    defineTypeNameAndDebug(ClausiusClapeyron, 0);
    addToRunTimeSelectionTable(singleComponentSatProp,ClausiusClapeyron, components);
}

Foam::ClausiusClapeyron::ClausiusClapeyron
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    singleComponentSatProp
    (
        mesh,
        dict
    ),
    pSat0_(modelDict().get<scalar>("pSat0")),
    TSat0_(modelDict().get<scalar>("TSat0")),
    L0_(modelDict().get<scalar>("L0")),
    R_(modelDict().get<scalar>("R")),
    C0_(C0())
{
}

Foam::scalar
Foam::ClausiusClapeyron::C0()
{
    // ln(p) = C0 - L/(R T)
    return L0_/(R_*TSat0_) + log(pSat0_);
}

Foam::scalar
Foam::ClausiusClapeyron::satT(scalar p) const
{
    // Invert ln(p) = C0 - L/(R T):
    // T = L / (R (C0 - ln(p))).
    // The previous implementation used ln(p)-C0 and therefore returned a
    // negative saturation temperature for physically valid inputs.
    return L0_/(R_ * (C0_ - Foam::log(p)));
}

Foam::scalar
Foam::ClausiusClapeyron::satP(Foam::scalar T) const
{
    return Foam::exp(C0_ - L0_/(R_*T));
}

Foam::scalar
Foam::ClausiusClapeyron::satL(Foam::scalar p) const
{
    return L0_;
}

// ************************************************************************* //
