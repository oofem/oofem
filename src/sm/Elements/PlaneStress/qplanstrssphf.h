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

#ifndef qplanstrssphf_h
#define qplanstrssphf_h

#include "../sm/Elements/phasefieldelement.h"
#include "../sm/Elements/PlaneStress/qplanstrss.h"

#define _IFT_QPlaneStressPhF2d_Name "QPlaneStressPhF2d"

namespace oofem {
/**
 * This class implements an Quadratic isoparametric eight-node quadrilateral plane-
 * stress phase field finite element. Each node has 3 degrees of freedom.
 * TODO in development
 */
class QPlaneStressPhF2d : public QPlaneStress2d, public PhaseFieldElement
{
protected:

public:
    QPlaneStressPhF2d(int n, Domain *d);
    virtual ~QPlaneStressPhF2d() { }

    NLStructuralElement *giveElement() { return this; } 

    virtual int computeNumberOfDofs() { return 24; }
    virtual void giveDofManDofIDMask( int inode, IntArray &answer ) const;
    virtual void giveDofManDofIDMask_u( IntArray &answer );
    virtual void giveDofManDofIDMask_d( IntArray &answer );

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_QPlaneStressPhF2d_Name; }
    virtual const char *giveClassName() const { return "QPlaneStressPhF2d"; }


    virtual void computeStiffnessMatrix( FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep ) 
    {
        PhaseFieldElement :: computeStiffnessMatrix( answer, rMode, tStep );
    }
    virtual void giveInternalForcesVector( FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0 ) 
    {
        PhaseFieldElement :: giveInternalForcesVector( answer, tStep, useUpdatedGpRecord );
    }
protected:

    
};
} // end namespace oofem
#endif // qplanstrss_h
