/* $Header: /home/cvs/bp/oofem/sm/src/truss1d.h,v 1.8 2003/04/06 14:08:32 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   *******************************
//   *** CLASS LUMPEDMASSELEMENT ***
//   *******************************


#ifndef lumpedmasselement_h
#define lumpedmasselement_h


#include "structuralelement.h"

namespace oofem {

class LumpedMassElement : public StructuralElement
{
    /*
     * This class implements a simple lumped mass element. Its purpose is to introduce 
     * an additional mass (mass components or rotary inertias) into a node.
     * The mass element is defined by a single node.
     * At present, mass is defined in the nodal coordinate system.
     * The same element can be used to add an additional stifness if needed (Not yet implemented).
     */

protected:
  FloatArray components; ///Mass and moments of inertia corresponding to nodal DOFs

public:
    LumpedMassElement(int, Domain *);                         // constructor
    ~LumpedMassElement()   { }                                // destructor

    void          computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void          computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }
    void          computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) 
    { answer.resize(0,0);}
    void          computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
    { answer.resize(0,0);}
    void                  computeNonForceLoadVector(FloatArray &answer, TimeStep *, ValueModeType mode)
    { answer.resize(0);}
    void computeForceLoadVector(FloatArray & answer, TimeStep *, ValueModeType)
    { answer.resize(0);}
    void giveInternalForcesVector(FloatArray &answer, TimeStep *, int useUpdatedGpRecord = 0)
    { answer.resize(0);}

    int            computeNumberOfDofs(EquationID ut);
    void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

    void                  updateInternalState(TimeStep *) {}
    void          updateYourself(TimeStep *tStep) {}
    int    checkConsistency();
#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    void          drawScalar(oofegGraphicContext &context);
#endif
    //
    // definition & identification
    //
    const char *giveClassName() const { return "LumpedMassElement"; }
    classType            giveClassID() const { return LumpedMassElementClass; }
    IRResultType initializeFrom(InputRecord *ir);
    Element_Geometry_Type giveGeometryType() const { return EGT_point; }

 protected:
    void  computeBmatrixAt(GaussPoint *, FloatMatrix &answer,
			   int lowerIndx = 1, int upperIndx = ALL_STRAINS)
    {}
    void  computeNmatrixAt(GaussPoint *, FloatMatrix &) {}
};

} // end namespace oofem
#endif // lumpedmasselement_h
