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

#ifndef Tr2Shell7XFEM_h
#define Tr2Shell7XFEM_h


#include "eleminterpmapperinterface.h"
#include "nodalaveragingrecoverymodel.h"
#include "layeredcrosssection.h"
#include "nlstructuralelement.h"
#include "shell7basexfem.h"

#define _IFT_Tr2Shell7XFEM_Name "tr2shell7xfem"

namespace oofem {

class FEI3dTrQuad;
class BoundaryLoad;

/**
 * This class represent a 7 parameter shell element. 
 * Each node has 7 degrees of freedom (displ. vec., director vec., inhomogeneous thickness strain ).
 * Add ref. to paper!
 * @author Jim Brouzoulis
 * @date 2012-11-01
 */

class Tr2Shell7XFEM : public Shell7BaseXFEM
{
protected:
    int numberOfGaussPoints;
    static FEI3dTrQuad interpolation;
    static bool __initialized;
    static IntArray ordering_phibar;
    static IntArray ordering_m;
    static IntArray ordering_gam;
    static IntArray ordering_all;
    static IntArray ordering_gr;
    static IntArray ordering_gr_edge;
    static bool initOrdering() {
        ordering_phibar.setValues(18, 1, 2, 3, 8, 9, 10, 15, 16, 17, 22, 23, 24, 29, 30 ,31, 36, 37, 38);
        ordering_m.setValues(18, 4, 5, 6, 11, 12, 13, 18, 19, 20, 25, 26, 27, 32, 33 ,34, 39, 40, 41);
        ordering_gam.setValues(6, 7, 14, 21, 28, 35, 42);
        ordering_all.setValues(42, 1, 2, 3, 8, 9, 10, 15, 16, 17, 22, 23, 24, 29, 30 ,31, 36, 37, 38,
                                   4, 5, 6, 11, 12, 13, 18, 19, 20, 25, 26, 27, 32, 33 ,34, 39, 40, 41,
                                   7, 14, 21, 28, 35, 42);
        ordering_gr.setValues(42, 1, 2, 3, 19, 20, 21, 37, 4, 5, 6, 22, 23, 24, 38, 7, 8, 9, 25, 26, 27, 39,
                              10, 11, 12, 28, 29, 30, 40, 13, 14, 15, 31, 32, 33, 41, 16, 17, 18,
                              34, 35, 36, 42);
        ordering_gr_edge.setValues(21, 1, 2, 3, 10, 11, 12, 19, 4, 5, 6, 13, 14, 15, 20, 7, 8, 9, 16, 17, 18, 21);
        return true;
    }

    virtual const IntArray &giveOrdering(SolutionField fieldType) const;

    //specific
    void giveSurfaceDofMapping(IntArray &answer, int iSurf) const;
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const;

    //virtual double computeVolumeAround(GaussPoint *gp);
    virtual double computeVolumeAroundLayer(GaussPoint *mastergp, int layer);
    virtual double computeAreaAround(GaussPoint *gp, double xi);

    virtual void computeGaussPoints();
    virtual void giveLocalNodeCoords(FloatArray &nodeLocalXiCoords, FloatArray &nodeLocalEtaCoords);

    virtual FEInterpolation *giveInterpolation();

    // VTK
    void vtkGiveUpdatedFictiousNodeCoords(FloatArray nodeCoords[15], int layer, TimeStep *tStep);

public:
    Tr2Shell7XFEM(int n, Domain *d);
    virtual ~Tr2Shell7XFEM() { }     // destructor -> declaring as virtual will make each subclass call their respective destr.
    // definition & identification
    //virtual int giveNumberOfDofs()           { return 42; }
    virtual int giveNumberOfEdgeDofs()       { return 21; }
    //virtual int giveNumberOfEdgeDofs()       { return 35; } ///@todo temporary! remove!
    virtual int giveNumberOfEdgeDofManagers(){ return 3;  }
    virtual const char *giveInputRecordName() const { return _IFT_Tr2Shell7XFEM_Name; }
    virtual const char *giveClassName() const { return "Tr2Shell7XFEM"; }
    virtual classType giveClassID() const { return Tr2Shell7XFEMClass; }
    //virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_2; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_Composite; }
    virtual integrationDomain giveIntegrationDomain() { return _Triangle; } // write new wedge-like type 'layeredWedge'
    virtual void giveCompositeExportData( IntArray &primaryVarsToExport, IntArray &cellVarsToExport,
                 std::vector<FloatArray> &nodeCoords, std::vector<IntArray> &cellNodes, IntArray &cellTypes, 
                 std::vector<FloatArray> &primaryVars, std::vector<FloatArray> &cellVars, TimeStep *tStep ){};

    //virtual void giveCompositeExportData( IntArray &primaryVarsToExport, IntArray &cellVarsToExport, TimeStep *tStep );
};


} // end namespace oofem
#endif 
