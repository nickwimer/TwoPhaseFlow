/*---------------------------------------------------------------------------*\
    Build wrapper for fixedSiteDetachedBubbleBirth.

    The implementation file defines its members with fully-qualified Foam::
    names and carries one terminal closing brace. Compile it inside an explicit
    C++ linkage block so that terminal brace closes this block without changing
    implementation logic.

    The temporary reduce wrapper below reports only the global quantities used
    by buildSiteStencils()'s startup validation. This keeps the diagnostic
    narrow while we qualify the new coarse model in parallel.
\*---------------------------------------------------------------------------*/

#include "fixedSiteDetachedBubbleBirth.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvPatch.H"
#include "fvScalarMatrix.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "surfaceFields.H"

#include <cmath>
#include <cstring>

template<class Type>
inline void fixedSiteStencilReduceDiagnostic
(
    const char* name,
    const Type& value
)
{
    if
    (
        Foam::Pstream::master()
     &&
        (
            std::strcmp(name, "localCaptureArea") == 0
         || std::strcmp(name, "localVolume") == 0
         || std::strcmp(name, "localMass") == 0
         || std::strcmp(name, "localLatent") == 0
         || std::strcmp(name, "localCoreCells") == 0
         || std::strcmp(name, "globalFaceCount") == 0
         || std::strcmp(name, "globalSphereCount") == 0
         || std::strcmp(name, "globalShieldCount") == 0
        )
    )
    {
        Foam::Info<< "FIXED_SITE_STENCIL_DIAGNOSTIC "
            << name << '=' << value << Foam::endl;
    }
}

#define reduce(value, operation)                                             \
    (Foam::reduce((value), (operation)),                                     \
     fixedSiteStencilReduceDiagnostic(#value, (value)))

extern "C++"
{
#include "fixedSiteDetachedBubbleBirth.C"
