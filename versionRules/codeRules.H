/*---------------------------------------------------------------------------* \
License
    This file is part of libWallModelledLES.

    libWallModelledLES is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    libWallModelledLES is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with libWallModelledLES. 
    If not, see <http://www.gnu.org/licenses/>.
 
Class


@brief
    Based on the OpenFOAM version, generates macros that can be used to make
    the code support multiple versions via #ifdef statements.
    The file is adopted with minor changes from swak4foam, see list of
    contributors below.

SourceFiles

Contributors/Copyright:
    2012-2018 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2015-2016 Bruno Santos <wyldckat@gmail.com>
    2016-2018 Mark Olesen <Mark.Olesen@esi-group.com>
    2019-2023 Timofey Mukha

\*---------------------------------------------------------------------------*/

#ifndef WmlesMacroHeader_H
#define WmlesMacroHeader_H

#include "foamVersion4wmles.H"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
#define VERSION_NR(major,minor,patch) (     \
    10000 * major + \
    100   * minor + \
    patch )

#define VERSION_NR2(major,minor) (     \
    10000 * major + \
    100   * minor + \
    99 )

#define FOAM_VERSION4WMLES VERSION_NR(FOAM_VERSION4WMLES_MAJOR,FOAM_VERSION4WMLES_MINOR,FOAM_VERSION4WMLES_PATCH_NUM)

// OpenFOAM-dev version. Has no number. Is assumed to be 2999.x
#if FOAM_VERSION4WMLES_MAJOR >= 2999
#define FOAM_OF_DEV_VERSION
#endif

#endif  // DOXYGEN_SHOULD_SKIP_THIS

#if (OPENFOAM_COM >= 1606) \
    || (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(5,0,0))
#define FOAM_HAS_TYPE_HEADER_OK
#if (OPENFOAM_COM >= 1712)
#define FOAM_TYPE_HEADER_OK_HAS_EXTRA_PARAMTERS
#endif
#endif

#if (OPENFOAM_COM >= 1606) \
    || (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(5,0,0))
#define FOAM_HAS_LOCAL_FILEPATH
#endif

#if (OPENFOAM_COM >= 1712) \
    || (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(5,0,0))
#define FOAM_LOCAL_FILEPATH_NEEDS_TYPE_PARAMETER
#endif

#if (OPENFOAM_COM >= 1712) \
    || (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(5,0,0))
#define FOAM_READSTREAM_METHOD_NEEDS_BOOL_PARAMETER
#endif

#if (OPENFOAM_COM >= 1606) \
    || (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>VERSION_NR2(3,0))
#define FOAM_COMPONENT_NAMES_ARE_CONST_PTR
#endif

#if (OPENFOAM_COM >= 1606) \
    || (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>VERSION_NR2(3,0))
#define FOAM_NEW_TMP_RULES
#endif

#if (OPENFOAM_COM >= 1612) \
    || (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>VERSION_NR2(3,0))
#define FOAM_PATCHFIELDTYPE_IN_GEOFIELD_IS_NOW_PATCH
#endif

#if (OPENFOAM_COM >= 1612) \
    || (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>VERSION_NR2(3,0))
// TODO: currently only used to comment the instantiations out
#define FOAM_CODED_FUNCTION_OBJECT_NEW_INTERFACE
#endif

#if (OPENFOAM_COM >= 1612) \
    || (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>VERSION_NR2(3,0))
#define FOAM_NEW_GEOMFIELD_RULES
#endif

#if (OPENFOAM_COM >= 1612) \
    || (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(5,0,0))
#define FOAM_MAKE_TEMPLATE_PATCHTYPE_FIELD_USES_PRIMITIVES
#endif

#if (OPENFOAM_COM >= 1706) \
    || (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(5,0,0))
#define FOAM_PSTREAM_COMMSTYPE_IS_ENUMCLASS
#endif

#if (OPENFOAM_COM >= 1712) \
    || (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(5,0,0))
#define FOAM_REGIOOBJECT_WRITEOBJECT_WITH_VALID
#endif

#if (OPENFOAM_COM >= 1803)
#define FOAM_LIST_HAS_MOVABLE_CONSTRUCT
#endif

#if (OPENFOAM_COM >= 1803)
#define FOAM_HASH_PTR_LIST_ACCEPTS_NO_RAW_POINTERS
#endif

#if (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(6,0,0))
#define FOAM_TIME_HAS_SEPARATE_SETDELTAT_NOADJUST
#endif

#if (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(6,0,0))
#define FOAM_TREEBOUNDBOX_DOES_NOT_ACCEPT_RNG
#endif

#if FOAM_VERSION4WMLES >= VERSION_NR(3,0,0) && !defined(FOAM_VERSION4WMLES_IS_EXTEND)
#define FOAM_NEW_TURBULENCE_STRUCTURE
#endif

#if (defined(FOAM_VERSION4WMLES_IS_COM))
#define FOAM_WALLDIST_CONSTRUCTOR_ACCEPTS_METHOD
#endif

#if (defined(FOAM_VERSION4WMLES_IS_COM) && OPENFOAM_COM >=1806)
/*! \def FOAM_AUTOPTR_HAS_CLONE_METHOD
    \brief Whether the autoPtr<> class has a clone() mehtod.
*/
#define FOAM_AUTOPTR_HAS_CLONE_METHOD
#endif

#if (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(7,0,0))
#define FOAM_VOLUMETYPE_NOT_CAPITAL
#endif

#if (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(7,0,0))
#define FOAM_NEW_WRITEENTRY
#endif

#if (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(8,0,0))
#define FOAM_DOTORG_TURBULENCE
#endif

#if (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(9,0,0))
#define FOAM_FVPATCHFIELD_NO_COPY
#endif

#if (defined(FOAM_VERSION4WMLES_IS_ORG))
#define FOAM_DICTIONARY_NO_GET
#endif

#if (defined(FOAM_VERSION4WMLES_IS_ORG) && FOAM_VERSION4WMLES>=VERSION_NR(8,0,0))
#define FOAM_DICTIONARY_HAS_LOOKUP
#endif

#if (defined(FOAM_VERSION4SWAK_IS_ORG) && FOAM_VERSION4SWAK >= VERSION_NR(9, 0, 0))
#define FOAM_FLUID_MOMENTUM_TRANSPORT_IS_NOW_DYNAMIC
#endif

#if (defined(FOAM_VERSION4SWAK_IS_ORG) && FOAM_VERSION4SWAK >= VERSION_NR(10, 0, 0))
#define FOAM_FLUID_MOMENTUM_TRANSPORT_IS_NO_LONGER_DYNAMIC
#endif

#endif

// ************************************************************************* //
