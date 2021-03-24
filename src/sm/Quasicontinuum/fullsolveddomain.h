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

#ifndef qc_fullsolveddomain_h
#define qc_fullsolveddomain_h

#include "floatarray.h"
#include "element.h"
#include "node.h"

#define _IFT_FullSolvedDomain_nodes "fsd_n"
#define _IFT_FullSolvedDomain_elements "fsd_e"
#define _IFT_FullSolvedDomain_radius "fsd_r"
#define _IFT_FullSolvedDomain_box "fsd_b"


//#define _IFT_QCFullsolveddomain_Name "QCFullsolveddomain"

namespace oofem {

/**
 * Information about fullsolved domain in CQ simulation.
 * 
 */
class QCFullsolveddomain 

{
protected:
    FloatArray FullSolvedDomainNodes;
    FloatArray FullSolvedDomainElements;
    FloatArray FullSolvedDomainRadius;
    FloatArray FullSolvedDomainBox;

public:
    QCFullsolveddomain();
    virtual ~QCFullsolveddomain();

    virtual void initializeFrom(InputRecord &ir);

    virtual void updateYourself();

    virtual bool isNodeInside(Node *n);


    // identification
    //virtual const char *giveInputRecordName() const { return _IFT_QCFullsolveddomain_Name; }
    virtual const char *giveClassName() const { return "QCFullsolveddomain"; }



};
} // end namespace oofem
#endif // qc_fullsolveddomain_h
