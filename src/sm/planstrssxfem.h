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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef planstrssxfem_h
#define planstrssxfem_h

#include "planstrss.h"
#include "xfemelementinterface.h"

namespace oofem {
/**
 * Temporary class for testing
 * in the usual case instead of PlaneStress2dXfem
 * there will be the standard PlaneStress2d
 */
class PlaneStress2dXfem : public PlaneStress2d, public XfemElementInterface
{
public:
    /// Constructor
    PlaneStress2dXfem(int n, Domain *d) : PlaneStress2d(n, d), XfemElementInterface(this) { }
    /// Destructor
    virtual ~PlaneStress2dXfem() { };

    virtual Interface *giveInterface(InterfaceType it);

    /// Computes the enriched part of the location array
    virtual void giveLocationArray(IntArray & locationArray, EquationID, const UnknownNumberingScheme & s) const;
    virtual const char *giveClassName() const { return "PlaneStress2dXfem"; }
    virtual classType giveClassID() const { return PlaneStress2dXfemClass; }
    virtual int computeNumberOfDofs(EquationID ut);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                          int lowerIndx = 1, int upperIndx = ALL_STRAINS);
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray & answer) const;
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *, TimeStep *tStep);
    virtual void computeVectorOf(EquationID type, ValueModeType u, TimeStep *stepN, FloatArray &answer);
    virtual void computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    //void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void drawScalar(oofegGraphicContext &context);
    //virtual void drawSpecial(oofegGraphicContext &);
    //void drawInternalState(oofegGraphicContext&);
#endif
};
} // end namespace oofem
#endif
