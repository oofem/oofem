/* $Header: /home/cvs/bp/oofem/oofemlib/src/element.h,v 1.27 2003/04/06 14:08:24 bp Exp $ */
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

/* Author: L. Svoboda */

//   ************************
//   *** CLASS TR01_SHELL ***
//   ************************
//   1.6.2010

#ifndef tr_shell01_h
#define tr_shell01_h

#include "structuralelement.h"

#include "cct3d.h"
#include "trplanrot3d.h"


namespace oofem {
class TR_SHELL01 : public StructuralElement
{
    /*
     * This class implements an triangular three-node shell finite element, composed of
     * cct3d and trplanrot3d elements.
     * Each node has 6 degrees of freedom.
     *
     */

protected:
    /// Pointer to plate element
    CCTPlate3d *plate;
    /// pointer to membrane (plane stress) element
    TrPlaneStrRot3d *membrane;

public:
    /// Constructor
    TR_SHELL01(int, Domain *);
    /// Destructor
    ~TR_SHELL01() { delete plate;
                    delete membrane; }

    virtual int computeNumberOfDofs(EquationID ut) { return 18; }
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
    { plate->giveDofManDofIDMask(inode, ut, answer); }

    //
    //// definition & identification
    //
    const char *giveClassName() const { return "TR_SHELL01"; }
    classType   giveClassID()   const { return TR_SHELL01Class; }
    IRResultType initializeFrom(InputRecord *ir);

    void giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep);
    void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep);

    void updateYourself(TimeStep *tStep);
    void updateInternalState(TimeStep *stepN);
    void printOutputAt(FILE *file, TimeStep *tStep);

    //
    // io routines
    //
#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    void          drawDeformedGeometry(oofegGraphicContext &, UnknownType type);
    //virtual void  drawScalar(oofegGraphicContext &context);
    //void          drawInternalState (oofegGraphicContext&);
#endif

protected:
    void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS)
    { _error("TR_SHELL01 :: computeBmatrixAt: calling of this function is not allowed"); }
    void computeNmatrixAt(GaussPoint *, FloatMatrix &)
    { _error("TR_SHELL01 :: computeNmatrixAt: calling of this function is not allowed"); }


    /// casem smazat
protected:
    integrationDomain giveIntegrationDomain()
    { return _Triangle;}
    void computeGaussPoints()
    { _error("TR_SHELL01 :: computeGaussPoints: calling of this function is not allowed"); }
    void computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
    { _error("TR_SHELL01 :: computeStressVector: calling of this function is not allowed"); }
    void computePrescribedStrainLocalLoadVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
    { _error("TR_SHELL01 :: ...: calling of this function is not allowed"); }
    void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode)
    { _error("TR_SHELL01 :: ...: calling of this function is not allowed"); }
    void computeForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
    { _error("TR_SHELL01 :: ...: calling of this function is not allowed"); }
    void computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
    { _error("TR_SHELL01 :: ...: calling of this function is not allowed"); }
public:
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
    { _error("TR_SHELL01 :: ...: calling of this function is not allowed"); }
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { _error("TR_SHELL01 :: ...: calling of this function is not allowed"); }
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
    { _error("TR_SHELL01 :: ...: calling of this function is not allowed"); }
};
} // end namespace oofem
#endif
