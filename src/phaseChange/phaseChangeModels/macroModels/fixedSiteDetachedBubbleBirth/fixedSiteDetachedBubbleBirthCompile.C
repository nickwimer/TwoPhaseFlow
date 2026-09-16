/*---------------------------------------------------------------------------*\
    Build wrapper for fixedSiteDetachedBubbleBirth.

    The implementation file defines its members with fully-qualified Foam::
    names and carries one terminal closing brace.  Compile it inside an explicit
    C++ linkage block so that terminal brace closes this block without changing
    any implementation logic.  This wrapper can be removed once the source is
    normalized to remove that terminal brace directly.
\*---------------------------------------------------------------------------*/

extern "C++"
{
#include "fixedSiteDetachedBubbleBirth.C"
