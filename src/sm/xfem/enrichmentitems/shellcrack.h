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
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef SHELLCRACK_H_
#define SHELLCRACK_H_

//#include "xfem/enrichmentitem.h"
#include "xfem/enrichmentitems/crack.h"

#define _IFT_ShellCrack_Name "shellcrack"
#define _IFT_ShellCrack_xiBottom "xibottom"
#define _IFT_ShellCrack_xiTop "xitop"

/**
 * Crack.
 * @author Jim Brouzoulis
 * @date July 29, 2014
 */
namespace oofem {
class XfemManager;
class Domain;
class InputRecord;
class GaussPoint;
class GnuplotExportModule;

class OOFEM_EXPORT ShellCrack : public Crack
{
public:
    ShellCrack(int n, XfemManager *xm, Domain *aDomain);

    virtual const char *giveClassName() const { return "ShellCrack"; }
    virtual const char *giveInputRecordName() const { return _IFT_ShellCrack_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    double xiBottom;
    double xiTop;
};
} // end namespace oofem

#endif /* SHELLCRACK_H_ */
