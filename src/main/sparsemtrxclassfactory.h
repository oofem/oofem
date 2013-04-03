/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "sparsemtrx.h"
#include "skyline.h"
#include "skylineu.h"
#include "spoolessparsemtrx.h"
#include "petscsparsemtrx.h"
#include "dss.h"
#ifdef __IML_MODULE
 #include "iml/compcol.h"
 #include "iml/dyncompcol.h"
 #include "iml/symcompcol.h"
 #include "iml/dyncomprow.h"
#endif

REGISTER_CLASS(Skyline, SMT_Skyline)
REGISTER_CLASS(SkylineUnsym, SMT_SkylineU)

#ifdef __IML_MODULE
REGISTER_CLASS(CompCol, SMT_CompCol)
REGISTER_CLASS(DynCompCol, SMT_DynCompCol)
REGISTER_CLASS(SymCompCol, SMT_SymCompCol)
REGISTER_CLASS(DynCompRow, SMT_DynCompRow)
#endif
#ifdef __SPOOLES_MODULE
REGISTER_CLASS(SpoolesSparseMtrx, SMT_SpoolesMtrx)
#endif
#ifdef __PETSC_MODULE
REGISTER_CLASS(PetscSparseMtrx, SMT_PetscMtrx)
#endif
#ifdef __DSS_MODULE
REGISTER_CLASS(DSSMatrixLDL, SMT_DSS_sym_LDL)
REGISTER_CLASS(DSSMatrixLL, SMT_DSS_sym_LL)
REGISTER_CLASS(DSSMatrixLU, SMT_DSS_unsym_LU)
#endif

