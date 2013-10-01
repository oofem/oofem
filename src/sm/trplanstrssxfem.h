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

#ifndef TRPLANSTRSSXFEM_H_
#define TRPLANSTRSSXFEM_H_

#include "trplanstrss.h"
#include "xfemelementinterface.h"
#include "vtkxmlexportmodule.h"


#define _IFT_TrPlaneStress2dXFEM_Name "trplanestress2dxfem"

namespace oofem {

/*
 *		trplanstrssxfem.h
 *
 *
 *      Class: TrPlaneStress2dXFEM
 *      Description: 3-node triangle with XFEM kinematics
 * 		@author Erik Svenning
 *
 */

class TrPlaneStress2dXFEM : public TrPlaneStress2d, public XfemElementInterface, public VTKXMLExportModuleElementInterface
{
public:

	TrPlaneStress2dXFEM(int n, Domain *d) : TrPlaneStress2d(n, d), XfemElementInterface(this), VTKXMLExportModuleElementInterface() { }

	virtual ~TrPlaneStress2dXFEM();


	virtual int checkConsistency();

#if 0
	// Overloaded functions from XfemElementInterface
    /// Partitions the element into patches by a triangulation.
    virtual void XfemElementInterface_partitionElement(AList< Triangle > *answer, std :: vector< FloatArray > &together);
    /// Updates integration rule based on the triangulation.
    virtual void XfemElementInterface_updateIntegrationRule();
    /// Helpful routine to put the nodes for triangulation together, should be in protected members probably.
    virtual void XfemElementInterface_prepareNodesForDelaunay(std::vector< std::vector< FloatArray > > &oPointPartitions);
#endif

    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }

    virtual Interface *giveInterface(InterfaceType it);

    virtual const char *giveInputRecordName() const { return _IFT_TrPlaneStress2dXFEM_Name; }
    virtual const char *giveClassName() const { return "TrPlaneStress2dXFEM"; }
    virtual classType giveClassID() const { return TrPlaneStress2dXFEMClass; }

    virtual int computeNumberOfDofs();
    virtual void computeGaussPoints();
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                          int lowerIndx = 1, int upperIndx = ALL_STRAINS);
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray & answer) const;
//    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *, TimeStep *tStep);
//    virtual void computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
//    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
//    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);

//    virtual Element_Geometry_Type giveGeometryType() const;

#ifdef __OOFEG
    // TODO: Implement OOFEG functions
    void drawRawGeometry(oofegGraphicContext &);
    //void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void drawScalar(oofegGraphicContext &context);
    //virtual void drawSpecial(oofegGraphicContext &);
    //void drawInternalState(oofegGraphicContext&);
#endif

};

} /* namespace oofem */
#endif /* TRPLANSTRSSXFEM_H_ */
