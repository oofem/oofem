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


#include "mixedgradientpressureneumann.h"
#include "dofiditem.h"
#include "dofmanager.h"
#include "dof.h"
#include "valuemodetype.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "engngm.h"
#include "node.h"
#include "element.h"
#include "integrationrule.h"
#include "gaussintegrationrule.h"
#include "gausspnt.h"
#include "masterdof.h"
#include "usrdefsub.h" // For sparse matrix creation.
#include "sparsemtrxtype.h"
#include "../fm/line2boundaryelement.h"
#include "mathfem.h"

#include "sparsemtrx.h"
#include "sparselinsystemnm.h"

namespace oofem {

MixedGradientPressureNeumann :: MixedGradientPressureNeumann(int n, Domain *d) : MixedGradientPressureBC(n,d)
{
    int nsd = d->giveNumberOfSpatialDimensions();
    int components = nsd*nsd - 1;
    this->sigmaDev = new Node(1, d); // Node number lacks meaning here.
    for (int i = 0; i < components; i++) {
        // Just putting in X_i id-items since they don't matter.
        sigmaDev->appendDof(new MasterDof(i+1, sigmaDev, (DofIDItem)(X_1+i) ));
    }
}


MixedGradientPressureNeumann :: ~MixedGradientPressureNeumann()
{
    delete sigmaDev;
}


int MixedGradientPressureNeumann :: giveNumberOfInternalDofManagers()
{
    return 1;
}


DofManager *MixedGradientPressureNeumann :: giveInternalDofManager(int i)
{
    return this->sigmaDev;
}


void MixedGradientPressureNeumann :: fromDeviatoricBase2D(FloatArray &cartesian, FloatArray &deviatoric)
{
    cartesian.resize(3);
    cartesian.at(1) = deviatoric.at(1)/sqrt(2.0);
    cartesian.at(2) = -deviatoric.at(1)/sqrt(2.0);
    cartesian.at(3) = (deviatoric.at(2) + deviatoric.at(3))*0.5;
}


void MixedGradientPressureNeumann :: fromDeviatoricBase3D(FloatArray &cartesian, FloatArray &deviatoric)
{
    cartesian.resize(6);
    cartesian.at(1) = 2.0 * deviatoric.at(1)/sqrt(6.0);
    cartesian.at(2) = (-deviatoric.at(1) + deviatoric.at(2))/sqrt(2.0);
    cartesian.at(3) = (-deviatoric.at(1) - deviatoric.at(2))/sqrt(2.0);
    //
    cartesian.at(4) = (deviatoric.at(3) + deviatoric.at(6))*0.5;
    cartesian.at(5) = (deviatoric.at(4) + deviatoric.at(7))*0.5;
    cartesian.at(6) = (deviatoric.at(5) + deviatoric.at(8))*0.5;
}


void MixedGradientPressureNeumann :: fromDeviatoricBase2D(FloatMatrix &cartesian, FloatMatrix &deviatoric)
{
    cartesian.resize(3,3);
    // E1 = [1/sqrt(2), 1/sqrt(2), 0,0]'; E2 = [0,0,1,0]'; E3 = [0,0,0,1]';
    // C = E1*(E1'*e11 + E2'*e12 + E3'*e13) + E2*(E1'*e21 + E2'*e22 + E3'*e23) + E3*(E1'*e31 + E2'*e32 + E3'*e33);
    
    cartesian.at(1,1) =   deviatoric.at(1,1) * 0.5;
    cartesian.at(2,2) =   deviatoric.at(1,1) * 0.5;
    cartesian.at(1,2) = - deviatoric.at(1,1) * 0.5;
    cartesian.at(2,1) = - deviatoric.at(1,1) * 0.5;
    //
    cartesian.at(1,3) =  (deviatoric.at(1,2) + deviatoric.at(1,3))/sqrt(8.0);
    cartesian.at(2,3) = -(deviatoric.at(1,2) + deviatoric.at(1,3))/sqrt(8.0);
    cartesian.at(3,1) =  (deviatoric.at(2,1) + deviatoric.at(3,1))/sqrt(8.0);
    cartesian.at(3,2) = -(deviatoric.at(2,1) + deviatoric.at(3,1))/sqrt(8.0);
    //
    cartesian.at(3,3) = (deviatoric.at(2,2) + deviatoric.at(2,3) + 
                         deviatoric.at(3,2) + deviatoric.at(3,3))*0.25;
}


void MixedGradientPressureNeumann :: fromDeviatoricBase3D(FloatMatrix &cartesian, FloatMatrix &deviatoric)
{
    cartesian.resize(6,6);
    /*
    syms 
    e11 e12 e13 e14 e15 e16 e17 e18 
    e21 e22 e23 e24 e25 e26 e27 e28 
    e31 e32 e33 e34 e35 e36 e37 e38 
    e41 e42 e43 e44 e45 e46 e47 e48 
    e51 e52 e53 e54 e55 e56 e57 e58 
    e61 e62 e63 e64 e65 e66 e67 e68 
    e71 e72 e73 e74 e75 e76 e77 e78 
    e81 e82 e83 e84 e85 e86 e87 e88;
    C = ...
    E1*(E1'*e11 + E2'*e12 + E3'*e13 + E4'*e14 + E5'*e15 + E6'*e16 + E7'*e17 + E8'*e18) + ...
    E2*(E1'*e21 + E2'*e22 + E3'*e23 + E4'*e24 + E5'*e25 + E6'*e26 + E7'*e27 + E8'*e28) + ...
    E3*(E1'*e31 + E2'*e32 + E3'*e33 + E4'*e34 + E5'*e35 + E6'*e36 + E7'*e37 + E8'*e38) + ...
    E4*(E1'*e41 + E2'*e42 + E3'*e43 + E4'*e44 + E5'*e45 + E6'*e46 + E7'*e47 + E8'*e48) + ...
    E5*(E1'*e51 + E2'*e52 + E3'*e53 + E4'*e54 + E5'*e55 + E6'*e56 + E7'*e57 + E8'*e58) + ...
    E6*(E1'*e61 + E2'*e62 + E3'*e63 + E4'*e64 + E5'*e65 + E6'*e66 + E7'*e67 + E8'*e68) + ...
    E7*(E1'*e71 + E2'*e72 + E3'*e73 + E4'*e74 + E5'*e75 + E6'*e76 + E7'*e77 + E8'*e78) + ...
    E8*(E1'*e81 + E2'*e82 + E3'*e83 + E4'*e84 + E5'*e85 + E6'*e86 + E7'*e87 + E8'*e88);
    */

    cartesian.at(1,1) = deviatoric.at(1,1)*2.0/3.0;
    cartesian.at(1,2) = - deviatoric.at(1,1)/3.0 + deviatoric.at(1,2)*sqrt(3.0);
    cartesian.at(2,1) = - deviatoric.at(1,1)/3.0 + deviatoric.at(1,2)*sqrt(3.0);
    cartesian.at(1,3) = - deviatoric.at(1,1)/3.0 - deviatoric.at(1,2)*sqrt(3.0);
    cartesian.at(3,1) = - deviatoric.at(1,1)/3.0 - deviatoric.at(1,2)*sqrt(3.0);

    cartesian.at(2,2) = deviatoric.at(1,1)/6.0 + deviatoric.at(2,2)/2.0 - deviatoric.at(1,2)/sqrt(12.0) - deviatoric.at(2,1)/sqrt(12.0);
    cartesian.at(2,3) = deviatoric.at(1,1)/6.0 - deviatoric.at(2,2)/2.0 + deviatoric.at(1,2)/sqrt(12.0) - deviatoric.at(2,1)/sqrt(12.0);
    cartesian.at(3,3) = deviatoric.at(1,1)/6.0 + deviatoric.at(2,2)/2.0 - deviatoric.at(1,2)/sqrt(12.0) + deviatoric.at(2,1)/sqrt(12.0);
    cartesian.at(3,2) = deviatoric.at(1,1)/6.0 - deviatoric.at(2,2)/2.0 + deviatoric.at(1,2)/sqrt(12.0) + deviatoric.at(2,1)/sqrt(12.0);
    // upper off diagonal part
    cartesian.at(1,4) = (deviatoric.at(1,3) + deviatoric.at(1,6))/sqrt(6.0);
    cartesian.at(1,5) = (deviatoric.at(1,4) + deviatoric.at(1,7))/sqrt(6.0);
    cartesian.at(1,6) = (deviatoric.at(1,5) + deviatoric.at(1,8))/sqrt(6.0);
    cartesian.at(2,4) = deviatoric.at(2,3)/sqrt(8.0) + deviatoric.at(2,6)/sqrt(8.0) - deviatoric.at(1,3)/sqrt(24.0) - deviatoric.at(1,6)/sqrt(24.0);
    cartesian.at(2,5) = deviatoric.at(2,4)/sqrt(8.0) + deviatoric.at(2,7)/sqrt(8.0) - deviatoric.at(1,4)/sqrt(24.0) - deviatoric.at(1,7)/sqrt(24.0);
    cartesian.at(2,6) = deviatoric.at(2,5)/sqrt(8.0) + deviatoric.at(2,8)/sqrt(8.0) - deviatoric.at(1,5)/sqrt(24.0) - deviatoric.at(1,8)/sqrt(24.0);
    cartesian.at(3,4) = - deviatoric.at(2,3)/sqrt(8.0) - deviatoric.at(2,6)/sqrt(8.0) - deviatoric.at(1,3)/sqrt(24.0) - deviatoric.at(1,6)/sqrt(24.0);
    cartesian.at(3,5) = - deviatoric.at(2,4)/sqrt(8.0) - deviatoric.at(2,7)/sqrt(8.0) - deviatoric.at(1,4)/sqrt(24.0) - deviatoric.at(1,7)/sqrt(24.0);
    cartesian.at(3,6) = - deviatoric.at(2,5)/sqrt(8.0) - deviatoric.at(2,8)/sqrt(8.0) - deviatoric.at(1,5)/sqrt(24.0) - deviatoric.at(1,8)/sqrt(24.0);
    // lower off diagonal part
    cartesian.at(1,4) = (deviatoric.at(3,1) + deviatoric.at(6,1))/sqrt(6.0);
    cartesian.at(1,5) = (deviatoric.at(4,1) + deviatoric.at(7,1))/sqrt(6.0);
    cartesian.at(1,6) = (deviatoric.at(5,1) + deviatoric.at(8,1))/sqrt(6.0);
    cartesian.at(2,4) = deviatoric.at(3,2)/sqrt(8.0) + deviatoric.at(6,2)/sqrt(8.0) - deviatoric.at(3,1)/sqrt(24.0) - deviatoric.at(6,1)/sqrt(24.0);
    cartesian.at(2,5) = deviatoric.at(4,2)/sqrt(8.0) + deviatoric.at(7,2)/sqrt(8.0) - deviatoric.at(4,1)/sqrt(24.0) - deviatoric.at(7,1)/sqrt(24.0);
    cartesian.at(2,6) = deviatoric.at(5,2)/sqrt(8.0) + deviatoric.at(8,2)/sqrt(8.0) - deviatoric.at(5,1)/sqrt(24.0) - deviatoric.at(8,1)/sqrt(24.0);
    cartesian.at(3,4) = - deviatoric.at(3,2)/sqrt(8.0) - deviatoric.at(6,2)/sqrt(8.0) - deviatoric.at(3,1)/sqrt(24.0) - deviatoric.at(6,1)/sqrt(24.0);
    cartesian.at(3,5) = - deviatoric.at(4,2)/sqrt(8.0) - deviatoric.at(7,2)/sqrt(8.0) - deviatoric.at(4,1)/sqrt(24.0) - deviatoric.at(7,1)/sqrt(24.0);
    cartesian.at(3,6) = - deviatoric.at(5,2)/sqrt(8.0) - deviatoric.at(8,2)/sqrt(8.0) - deviatoric.at(5,1)/sqrt(24.0) - deviatoric.at(8,1)/sqrt(24.0);
    //
    cartesian.at(4,4) = (deviatoric.at(3,3) + deviatoric.at(3,6) + deviatoric.at(6,3) + deviatoric.at(6,6))*0.25;
    cartesian.at(4,5) = (deviatoric.at(3,4) + deviatoric.at(3,7) + deviatoric.at(6,4) + deviatoric.at(6,7))*0.25;
    cartesian.at(4,6) = (deviatoric.at(3,5) + deviatoric.at(3,8) + deviatoric.at(6,5) + deviatoric.at(6,8))*0.25;
    cartesian.at(5,4) = (deviatoric.at(4,3) + deviatoric.at(4,6) + deviatoric.at(7,3) + deviatoric.at(7,6))*0.25;
    cartesian.at(5,5) = (deviatoric.at(4,4) + deviatoric.at(4,7) + deviatoric.at(7,4) + deviatoric.at(7,7))*0.25;
    cartesian.at(5,6) = (deviatoric.at(4,5) + deviatoric.at(4,8) + deviatoric.at(7,5) + deviatoric.at(7,8))*0.25;
    cartesian.at(6,4) = (deviatoric.at(5,3) + deviatoric.at(5,6) + deviatoric.at(8,3) + deviatoric.at(8,6))*0.25;
    cartesian.at(6,5) = (deviatoric.at(5,4) + deviatoric.at(5,7) + deviatoric.at(8,4) + deviatoric.at(8,7))*0.25;
    cartesian.at(6,6) = (deviatoric.at(5,5) + deviatoric.at(5,8) + deviatoric.at(8,5) + deviatoric.at(8,8))*0.25;
}


void MixedGradientPressureNeumann :: setPrescribedDeviatoricGradientFromVoigt(const FloatArray &t)
{
    // Converts the Voigt vector to a deviatoric base dyads (which don't assume symmetry)
    int nsd = this->domain->giveNumberOfSpatialDimensions();
    if (nsd == 3) {
        this->devGradient.resize(8);
        this->devGradient.at(1) = (2.0*t.at(1) - t.at(2) - t.at(3))/sqrt(6.);
        this->devGradient.at(2) = (t.at(2) - t.at(3))/sqrt(2.);
        //
        this->devGradient.at(3) = 0.5*t.at(4);
        this->devGradient.at(4) = 0.5*t.at(5);
        this->devGradient.at(5) = 0.5*t.at(6);
        //
        this->devGradient.at(6) = 0.5*t.at(4);
        this->devGradient.at(7) = 0.5*t.at(5);
        this->devGradient.at(8) = 0.5*t.at(6);
        
        this->volGradient = t.at(1) + t.at(2) + t.at(3);
    } else if (nsd == 2) {
        this->devGradient.resize(3);
        this->devGradient.at(1) = (t.at(1) - t.at(2))/sqrt(2.);
        this->devGradient.at(2) = 0.5*t.at(3);
        this->devGradient.at(3) = 0.5*t.at(3);
        
        this->volGradient = t.at(1) + t.at(2);
    } else {
        this->devGradient.resize(0);
        this->volGradient = t.at(1);
    }
}


void MixedGradientPressureNeumann :: giveLocationArrays(AList<IntArray> &rows, AList<IntArray> &cols, EquationID eid, CharType type,
    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, Domain *domain)
{
    if (eid == EID_MomentumBalance_ConservationEquation)
        eid = EID_MomentumBalance;

    if (eid != EID_MomentumBalance || type != TangentStiffnessMatrix) {
        return;
    }

    IntArray loc_r, loc_c, sigma_loc_r, sigma_loc_c;

    // Fetch the columns/rows for the stress contributions;
    this->sigmaDev->giveCompleteLocationArray(sigma_loc_r, r_s);
    this->sigmaDev->giveCompleteLocationArray(sigma_loc_c, c_s);

    rows.growTo(this->boundaries.size()*2);
    cols.growTo(this->boundaries.size()*2);
    int i = 1;
    for (dynaList< std::pair<int,int> > :: iterator pos = this->boundaries.begin(); pos != this->boundaries.end(); ++pos) {
        Element *e = this->giveDomain()->giveElement( (*pos).first );
        int boundary = (*pos).second;

        e->giveBoundaryLocationArray(loc_r, boundary, eid, r_s);
        e->giveBoundaryLocationArray(loc_c, boundary, eid, c_s);
        // For most uses, loc_r == loc_c, and sigma_loc_r == sigma_loc_c.
        rows.put(i, new IntArray(loc_r));
        cols.put(i, new IntArray(sigma_loc_c));
        i++;
        // and the symmetric part (usually the transpose of above)
        rows.put(i, new IntArray(sigma_loc_r));
        cols.put(i, new IntArray(loc_c));
        i++;
    }
}


IntegrationRule *MixedGradientPressureNeumann :: CreateIntegrationRule(Element *e, int order)
{
    // The element should give us most/all of this information;
    GaussIntegrationRule *ir = new GaussIntegrationRule(1, e);
    MaterialMode matMode = e->giveMaterialMode();
    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();
    int npoints;
    integrationDomain id;
    if (nsd == 3) {
        npoints = 4; ///@todo I don't know how to determine this for surfaces
        id = _Triangle; ///@todo We need to obtain this from the element itself.
    } else if (nsd == 2) {
        npoints = (order + 1 + 1)/2; // extra +1 for rounding up; npoints gives exact integration for order = npoints*2 - 1
        id = _Line;
    } else {
        npoints = 1;
        id = _Point;
    }
    ir->setUpIntegrationPoints(id, npoints, matMode);
    return ir;
}


void MixedGradientPressureNeumann :: integrateVolTangent(FloatArray &answer, Element *e, int boundary)
{
    FloatArray normal, n, contrib;
    FloatMatrix nMatrix;
    IntArray boundaryNodes;
    
    FEInterpolation *interp = e->giveInterpolation(); // Geometry interpolation
    // Some assumptions here, either velocity or displacement unknowns. Perhaps its better to just as for a certain equation id?
    // Code will be written using v for velocity, but it could represent displacements.
    FEInterpolation *interpUnknown = e->giveInterpolation(V_u);
    if (interpUnknown) {
        interpUnknown = e->giveInterpolation(D_u);
    }
    
    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();
    int order = interp->giveInterpolationOrder() + interpUnknown->giveInterpolationOrder();
    IntegrationRule *ir = this->CreateIntegrationRule(e, order);

    answer.resize(0);
    for (int i = 0; i < ir->getNumberOfIntegrationPoints(); i++) {
        GaussPoint *gp = ir->getIntegrationPoint(i);
        FloatArray &lcoords = *gp->giveCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        // Evaluate the normal;
        double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
        // Evaluate the velocity/displacement coefficients
        interpUnknown->boundaryEvalN(n, lcoords, cellgeo);
        nMatrix.beNMatrixOf(n, nsd);

        contrib.beTProductOf(nMatrix, normal);
        
        answer.add(detJ*gp->giveWeight(), contrib);
    }
}


void MixedGradientPressureNeumann :: integrateDevTangent(FloatMatrix &answer, Element *e, int boundary)
{
    FloatArray normal, n;
    FloatMatrix nMatrix, E_n;
    FloatMatrix contrib;
    
    FEInterpolation *interp = e->giveInterpolation(); // Geometry interpolation
    // Some assumptions here, either velocity or displacement unknowns. Perhaps its better to just as for a certain equation id?
    // Code will be written using v for velocity, but it could represent displacements.
    FEInterpolation *interpUnknown = e->giveInterpolation(V_u);
    if (interpUnknown) {
        interpUnknown = e->giveInterpolation(D_u);
    }
    
    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();
    int order = interp->giveInterpolationOrder() + interpUnknown->giveInterpolationOrder();
    IntegrationRule *ir = this->CreateIntegrationRule(e, order);
    
    answer.resize(0,0);
    for (int i = 0; i < ir->getNumberOfIntegrationPoints(); i++) {
        GaussPoint *gp = ir->getIntegrationPoint(i);
        FloatArray &lcoords = *gp->giveCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        // Evaluate the normal;
        double detJ = interp->boundaryEvalNormal(normal, boundary, lcoords, cellgeo);
        // Evaluate the velocity/displacement coefficients
        interpUnknown->boundaryEvalN(n, lcoords, cellgeo);
        nMatrix.beNMatrixOf(n, nsd);

        // Formulating like this to avoid third order tensors, which is hard to express in linear algebra.
        // v (x) n : E_i = v . ( E_i . n ) = v . E_n 
        if (nsd == 3) {
            E_n.resize(8,3);
            E_n.at(1,1) = 2.0*normal.at(1)/sqrt(6.0);
            E_n.at(1,2) = -normal.at(2)/sqrt(6.0);
            E_n.at(1,3) = -normal.at(3)/sqrt(6.0);
            
            E_n.at(2,1) = 0.;
            E_n.at(2,2) = normal.at(2)/sqrt(2.0);
            E_n.at(2,3) = -normal.at(3)/sqrt(2.0);
            
            E_n.at(3,1) = normal.at(2);
            E_n.at(3,2) = 0.;
            E_n.at(3,3) = 0.;
            
            E_n.at(4,1) = normal.at(3);
            E_n.at(4,2) = 0.;
            E_n.at(4,3) = 0.;
            
            E_n.at(5,1) = 0.;
            E_n.at(5,2) = normal.at(3);
            E_n.at(5,3) = 0.;

            E_n.at(6,1) = 0.;
            E_n.at(6,2) = normal.at(1);
            E_n.at(6,3) = 0.;

            E_n.at(7,1) = 0.;
            E_n.at(7,2) = 0.;
            E_n.at(7,3) = normal.at(1);

            E_n.at(8,1) = 0.;
            E_n.at(8,2) = 0.;
            E_n.at(8,3) = normal.at(2);

        } else if (nsd == 2) {
            E_n.resize(3,2);
            E_n.at(1,1) = normal.at(1)/sqrt(2.0);
            E_n.at(1,2) = -normal.at(2)/sqrt(2.0);
            
            E_n.at(2,1) = normal.at(2);
            E_n.at(2,2) = 0.;
            
            E_n.at(3,1) = 0.;
            E_n.at(3,2) = normal.at(1);
        } else {
            E_n.resize(0,0);
        }
        
        contrib.beProductOf(E_n, nMatrix);
        
        answer.add(detJ*gp->giveWeight(), contrib);
    }
}


double MixedGradientPressureNeumann :: assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                    CharType type, ValueModeType mode, const UnknownNumberingScheme &s, Domain *domain)
{
    // Boundary condition only acts on the momentumbalance part.
    if (eid == EID_MomentumBalance_ConservationEquation)
        eid = EID_MomentumBalance;

    if (eid != EID_MomentumBalance)
        return 0.0;

    double norm;
    IntArray loc, sigma_loc;  // For the velocities and stress respectively
    this->sigmaDev->giveCompleteLocationArray(sigma_loc, s);

    if (type == ExternalForcesVector) {
        // The external forces have two contributions. On the additional equations for sigmaDev, the load is simple the deviatoric gradient.
        double rve_size = this->domainSize();
        FloatArray devLoad;
        devLoad.beScaled(-rve_size, this->devGradient);
        answer.assemble(devLoad, sigma_loc);
        norm = devLoad.computeSquaredNorm(); ///@todo Different units here, have to consider something better..

        // The second contribution is on the momentumbalance equation; - int delta_v . n dA * p
        FloatArray fe;
        for (dynaList< std::pair<int,int> > :: iterator pos = this->boundaries.begin(); pos != this->boundaries.end(); ++pos) {
            Element *e = this->giveDomain()->giveElement( (*pos).first );
            int boundary = (*pos).second;
            
            e->giveBoundaryLocationArray(loc, boundary, eid, s);
            this->integrateVolTangent(fe, e, boundary);
            fe.times(-this->pressure);
            answer.assemble(fe, loc);
        }
        return norm;
        
    } else if (type == InternalForcesVector) {
        FloatMatrix Ke;
        FloatArray fe_v, fe_s;
        FloatArray s_dev, e_v;
        
        // Fetch the current values of the stress;
        s_dev.resize(this->sigmaDev->giveNumberOfDofs());
        for ( int i = 1; i <= this->sigmaDev->giveNumberOfDofs(); i++ ) {
            s_dev.at(i) = this->sigmaDev->giveDof(i)->giveUnknown(eid, mode, tStep);
        }
        // Assemble: int delta_v (x) n dA : E_i s_i
        //           int v (x) n dA : E_i delta_s_i
        norm = 0.;
        for (dynaList< std::pair<int,int> > :: iterator pos = this->boundaries.begin(); pos != this->boundaries.end(); ++pos) {
            Element *e = this->giveDomain()->giveElement( (*pos).first );
            int boundary = (*pos).second;
            
            // Fetch the element information;
            e->giveBoundaryLocationArray(loc, boundary, eid, s);
            e->computeBoundaryVectorOf(boundary, eid, mode, tStep, e_v);
            this->integrateDevTangent(Ke, e, boundary);
            
            // We just use the tangent, less duplicated code (the addition of sigmaDev is linear).
            fe_v.beProductOf(Ke, e_v);
            fe_s.beTProductOf(Ke, s_dev);
            
            answer.assemble(fe_s, loc); // Contributions to delta_v equations
            answer.assemble(fe_v, sigma_loc); // Contribution to delta_s_i equations
            norm += fe_v.computeNorm() /*+ fe_s.computeNorm()*/; ///@todo What to do with the mixed unknowns here? 
        }
        return norm;
    }
    return 0.;
}


void MixedGradientPressureNeumann :: assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
    CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, Domain *domain)
{
    if (eid == EID_MomentumBalance_ConservationEquation)
        eid = EID_MomentumBalance;

    if (eid != EID_MomentumBalance)
        return;

    if (type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == StiffnessMatrix || type == ElasticStiffnessMatrix) {
        FloatMatrix Ke, KeT;
        IntArray loc_r, loc_c, sigma_loc_r, sigma_loc_c;

        // Fetch the columns/rows for the stress contributions;
        this->sigmaDev->giveCompleteLocationArray(sigma_loc_r, r_s);
        this->sigmaDev->giveCompleteLocationArray(sigma_loc_c, c_s);
        
        for (dynaList< std::pair<int,int> > :: iterator pos = this->boundaries.begin(); pos != this->boundaries.end(); ++pos) {
            Element *e = this->giveDomain()->giveElement( (*pos).first );
            int boundary = (*pos).second;
            
            // Fetch the element information;
            e->giveBoundaryLocationArray(loc_r, boundary, eid, r_s);
            e->giveBoundaryLocationArray(loc_c, boundary, eid, c_s);
            this->integrateDevTangent(Ke, e, boundary);
            KeT.beTranspositionOf(Ke);

            answer->assemble(sigma_loc_r, loc_c, Ke); // Contribution to delta_s_i equations
            answer->assemble(loc_r, sigma_loc_c, KeT); // Contributions to delta_v equations
        }
    } 
}


void MixedGradientPressureNeumann :: computeFields(FloatArray &sigmaDev, double &vol, EquationID eid, TimeStep *tStep)
{
    int nsd = this->giveDomain()->giveNumberOfSpatialDimensions();
    FloatArray sigmaDevBase;
    IntArray dofMask(0);
    // Fetch the current values of the stress in deviatoric base;
    sigmaDevBase.resize(this->sigmaDev->giveNumberOfDofs());
    for ( int i = 1; i <= this->sigmaDev->giveNumberOfDofs(); i++ ) {
        sigmaDevBase.at(i) = this->sigmaDev->giveDof(i)->giveUnknown(eid, VM_Total, tStep);
    }
    // Convert it back from deviatoric base:
    if (nsd == 3) {
        this->fromDeviatoricBase3D(sigmaDev, sigmaDevBase);
    } else if (nsd == 2) {
        this->fromDeviatoricBase2D(sigmaDev, sigmaDevBase);
    } else {
        sigmaDev.resize(0);
    }

    // Postprocessing; vol = int v . n dA
    FloatArray unknowns, fe;
    vol = 0.;
    for (dynaList< std::pair<int,int> > :: iterator pos = this->boundaries.begin(); pos != this->boundaries.end(); ++pos) {
        Element *e = this->giveDomain()->giveElement( (*pos).first );
        int boundary = (*pos).second;
        
        e->computeVectorOf(eid, VM_Total, tStep, unknowns);
        this->integrateVolTangent(fe, e, boundary);
        vol += fe.dotProduct(unknowns);
    }
    double rve_size = this->domainSize();
    vol /= rve_size;
    vol -= volGradient; // This is needed for consistency; We return the volumetric "residual" if a gradient with volumetric contribution is set.
}


void MixedGradientPressureNeumann :: computeTangents(
    FloatMatrix &Ed, FloatArray &Ep, FloatArray &Cd, double &Cp, EquationID eid, TimeStep *tStep)
{
    //double size = this->domainSize();
    // Fetch some information from the engineering model
    EngngModel *rve = this->giveDomain()->giveEngngModel();
    ///@todo Get this from engineering model
    SparseLinearSystemNM *solver = CreateUsrDefSparseLinSolver(ST_Petsc, 1, this->domain, this->domain->giveEngngModel());// = rve->giveLinearSolver();
    SparseMtrx *Kff;
    SparseMtrxType stype = SMT_PetscMtrx;// = rve->giveSparseMatrixType();
    EModelDefaultEquationNumbering fnum;
    double rve_size = this->domainSize();
    
    // Set up and assemble tangent FE-matrix which will make up the sensitivity analysis for the macroscopic material tangent.
    Kff = CreateUsrDefSparseMtrx(stype);
    if ( !Kff ) {
        OOFEM_ERROR2("MixedGradientPressureNeumann :: computeTangents - Couldn't create sparse matrix of type %d\n", stype);
    }
    Kff->buildInternalStructure( rve, 1, eid, fnum, fnum );
    rve->assemble(Kff, tStep, eid, StiffnessMatrix, fnum, fnum, this->domain );

    // Setup up indices and locations
    int neq = Kff->giveNumberOfRows();

    // Indices and such of internal dofs
    int ndev = this->sigmaDev->giveNumberOfDofs();

    // Matrices and arrays for sensitivities
    FloatMatrix ddev_pert(neq,ndev); // In fact, npeq should most likely equal ndev
    FloatArray p_pert(neq); // RHS for d_dev [d_dev11, d_dev22, d_dev12] in 2D

    FloatMatrix s_d(neq,ndev); // Sensitivity fields for d_dev
    FloatArray s_p(neq); // Sensitivity fields for p

    // Unit pertubations for d_dev
    ddev_pert.zero();
    for (int i = 1; i <= ndev; ++i) {
        int eqn = this->sigmaDev->giveDof(i)->giveEquationNumber(fnum);
        ddev_pert.at(eqn, i) = -1.0*rve_size;
    }
    
    // Unit pertubation for d_p
    p_pert.zero();
    FloatArray fe;
    IntArray loc;
    for (dynaList< std::pair<int,int> > :: iterator pos = this->boundaries.begin(); pos != this->boundaries.end(); ++pos) {
        Element *e = this->giveDomain()->giveElement( (*pos).first );
        int boundary = (*pos).second;

        e->giveBoundaryLocationArray(loc, boundary, eid, fnum);
        this->integrateVolTangent(fe, e, boundary);
        fe.times(-1.0); // here d_p = 1.0
        p_pert.assemble(fe, loc);
    }

    // Solve all sensitivities
    solver->solve(Kff,ddev_pert,s_d);
    solver->solve(Kff,&p_pert,&s_p);

    // Extract the stress response from the solutions
    FloatArray sigma_p(ndev);
    FloatMatrix sigma_d(ndev,ndev);
    for (int i = 1; i <= ndev; ++i) {
        int eqn = this->sigmaDev->giveDof(i)->giveEquationNumber(fnum);
        sigma_p.at(i) = s_p.at(eqn);
        for (int j = 1; j <= ndev; ++j) {
            sigma_d.at(i,j) = s_d.at(eqn,j);
        }
    }

    // Post-process the volumetric rate of deformations in the sensitivity fields;
    FloatArray e_d(ndev); e_d.zero();
    double e_p = 0.0;
    for (dynaList< std::pair<int,int> > :: iterator pos = this->boundaries.begin(); pos != this->boundaries.end(); ++pos) {
        Element *e = this->giveDomain()->giveElement( (*pos).first );
        int boundary = (*pos).second;
    
        this->integrateVolTangent(fe, e, boundary);
        e->giveBoundaryLocationArray(loc, boundary, eid, fnum);
        
        // Using "loc" to pick out the relevant contributions. This won't work at all if there are local coordinate systems in these nodes
        // or slave nodes etc. The goal is to compute the velocity from the sensitivity field, but we need to avoid going through the actual
        // engineering model. If this ever becomes an issue it needs to perform the same steps as Element::giveUnknownVector does.
        for (int i = 1; i <= fe.giveSize(); ++i) {
            if (loc.at(i) > 0)  {
                e_p += fe.at(i) * s_p.at(loc.at(i));
                for (int j = 1; j <= ndev; ++j) {
                    e_d.at(j) += fe.at(i) * s_d.at(loc.at(i), j);
                }
            }
        }
    }
    e_p /= rve_size;

    // Now we need to express the tangents in the normal cartesian coordinate system (as opposed to the deviatoric base we use during computations
    Cp = e_p; // Scalar components are of course the same
    double nsd = this->giveDomain()->giveNumberOfSpatialDimensions();
    if (nsd == 3) {
        this->fromDeviatoricBase3D(Cd, e_d);
        this->fromDeviatoricBase3D(Ep, sigma_p);
        this->fromDeviatoricBase3D(Ed, sigma_d);
        
    } else if (nsd == 2) {
        this->fromDeviatoricBase2D(Cd, e_d);
        this->fromDeviatoricBase2D(Ep, sigma_p);
        this->fromDeviatoricBase2D(Ed, sigma_d);

    } else { // For 1D case, there simply are no deviatoric components!
        Cd.resize(0);
        Ep.resize(0);
        Ed.beEmptyMtrx();
    }

    delete Kff;    
    delete solver; ///@todo Remove this when solver is taken from engngmodel
}


IRResultType MixedGradientPressureNeumann :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;

    ActiveBoundaryCondition :: initializeFrom(ir);

    FloatArray gradient;
    IR_GIVE_FIELD(ir, gradient, IFT_MixedGradientPressure_devGradient, "devgradient");
    this->setPrescribedDeviatoricGradientFromVoigt(gradient);
    IR_GIVE_FIELD(ir, this->pressure, IFT_MixedGradientPressure_pressure, "pressure");

    return IRRT_OK;
}


int MixedGradientPressureNeumann :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    GeneralBoundaryCondition :: giveInputRecordString(str, keyword);

    sprintf( buff, " devgradient %d ", this->devGradient.giveSize() );
    for ( int i = 1; i <= this->devGradient.giveSize(); i++ ) {
        sprintf( buff, " %e", this->devGradient.at(i) );
        str += buff;
    }

    return 1;
}


void MixedGradientPressureNeumann :: addElementSide(int elem, int side)
{
    this->boundaries.pushBack(std::pair<int,int>(elem,side));
}


void MixedGradientPressureNeumann :: addElement(int elem)
{
    this->addElementSide(elem, 0);
}


void MixedGradientPressureNeumann :: clearElements()
{
    this->boundaries.clear();
}


void MixedGradientPressureNeumann :: scale(double s)
{
    devGradient.times(s);
    pressure *= s;
}

} // end namespace oofem

