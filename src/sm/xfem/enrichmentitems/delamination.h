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

#ifndef DELAMINATION_H_
#define DELAMINATION_H_

#include "xfem/enrichmentitem.h"

///@name Input fields for Delamination
//@{
#define _IFT_Delamination_Name "delamination"
#define _IFT_Delamination_xiCoord "delaminationxicoord"
#define _IFT_Delamination_interfacenum "interfacenum"
#define _IFT_Delamination_csnum "csnum"
#define _IFT_Delamination_CohesiveZoneMaterial "czmaterial"
//#define _IFT_MultipleDelamination_Name "multipledelamination"
//@}

namespace oofem {
class Material;

/**
 * Delamination.
 * */
class OOFEM_EXPORT Delamination : public EnrichmentItem
{
protected:
    Material *mat;  // Material for cohesive zone model
    int interfaceNum;
    int crossSectionNum;
    int matNum;
    double delamXiCoord;    // defines at what local xi-coord the delamination is defined
public:
    Delamination(int n, XfemManager *xm, Domain *aDomain);

    virtual const char *giveClassName() const { return "Delamination"; }
    virtual const char *giveInputRecordName() const { return _IFT_Delamination_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void appendInputRecords(DynamicDataReader &oDR);

    double giveDelamXiCoord() { return delamXiCoord; };
    //virtual Material *giveMaterial() { return mat; }
    virtual void updateGeometry(FailureCriteriaStatus *fc, TimeStep *tStep);
};
} // end namespace oofem


#endif /* DELAMINATION_H_ */
