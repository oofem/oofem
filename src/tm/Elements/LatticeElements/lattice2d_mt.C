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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "tm/Elements/LatticeElements/lattice2d_mt.h"
#include "tm/Materials/transportmaterial.h"
#include "tm/Materials/LatticeMaterials/latticetransmat.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "load.h"
#include "classfactory.h"
#include "parametermanager.h"
#include "paramkey.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
REGISTER_Element(Lattice2d_mt);

ParamKey Lattice2d_mt::IPK_Lattice2d_mt_dim("dim");
ParamKey Lattice2d_mt::IPK_Lattice2d_mt_thickness("thick");
ParamKey Lattice2d_mt::IPK_Lattice2d_mt_width("width");
ParamKey Lattice2d_mt::IPK_Lattice2d_mt_gpcoords("gpcoords");
ParamKey Lattice2d_mt::IPK_Lattice2d_mt_crackwidth("crackwidth");
ParamKey Lattice2d_mt::IPK_Lattice2d_mt_couplingflag("couplingflag");
ParamKey Lattice2d_mt::IPK_Lattice2d_mt_couplingnumber("couplingnumber");

Lattice2d_mt :: Lattice2d_mt(int n, Domain *aDomain, ElementMode em) :
    LatticeTransportElement(n, aDomain, em)
{
    numberOfDofMans  = 2;
    dimension = 2.;
    numberOfGaussPoints = 1;
    couplingNumbers.resize(1);
    couplingNumbers.zero();
    crackWidths.resize(1);
    crackWidths.zero();
}

double Lattice2d_mt :: giveLength()
{
    if ( length == 0. ) {
        Node *nodeA = this->giveNode(1);
        Node *nodeB = this->giveNode(2);
        double dx   = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        double dy   = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
        length = sqrt(dx * dx + dy * dy);
    }

    return length;
}

double
Lattice2d_mt :: givePressure()
{
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    GaussPoint *gp = iRule->getIntegrationPoint(0);
    auto status = static_cast< LatticeTransportMaterialStatus * >( gp->giveMaterialStatus() );

    return status->givePressure();
}

double
Lattice2d_mt :: giveOldPressure()
{
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    GaussPoint *gp = iRule->getIntegrationPoint(0);
    auto status = static_cast< LatticeTransportMaterialStatus * >( gp->giveMaterialStatus() );

    return status->giveOldPressure();
}


double
Lattice2d_mt :: giveMass()
{
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    GaussPoint *gp = iRule->getIntegrationPoint(0);
    auto status = static_cast< LatticeTransportMaterialStatus * >( gp->giveMaterialStatus() );

    double mass = status->giveMass();
    //multiply with volume
    mass *= this->length * this->width / 2.;
    return mass;
}


void
Lattice2d_mt :: computeNSubMatrixAt(FloatMatrix &answer, const FloatArray &coords)
{
    double ksi = coords.at(1);
    double n1 = ( 1. - ksi ) * 0.5;
    double n2 = ( 1. + ksi ) * 0.5;

    answer.resize(1, 2);
    answer.zero();

    answer.at(1, 1) = n1;
    answer.at(1, 2) = n2;
}

void
Lattice2d_mt :: computeNmatrixAt(FloatMatrix &answer, const FloatArray &coords)
{
    this->computeNSubMatrixAt(answer, coords);
}


void
Lattice2d_mt :: computeGradientMatrixAt(FloatMatrix &answer, const FloatArray &lcoords)
{
    double l = this->giveLength();

    //    Assemble Gradient Matrix used to compute temperature gradient
    answer.resize(1, 2);
    answer.zero();
    answer.at(1, 1) = -1.;
    answer.at(1, 2) = 1.;

    answer.times(1. / l);
}

void
Lattice2d_mt :: updateInternalState(TimeStep *tStep)
{
    FloatArray f, r;
    FloatMatrix n;
    TransportMaterial *mat = static_cast< TransportMaterial * >( this->giveMaterial() );

    // force updating ip values
    for ( auto &iRule: integrationRulesArray ) {
        for ( auto &gp: * iRule ) {
            this->computeNmatrixAt(n, gp->giveNaturalCoordinates() );
            this->computeVectorOf({ P_f }, VM_Total, tStep, r);
            f.beProductOf(n, r);
            mat->updateInternalState(f, gp, tStep);
        }
    }
}


void
Lattice2d_mt :: computeGaussPoints()
{
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ] = std :: make_unique< GaussIntegrationRule >(1, this, 1, 2);
    integrationRulesArray [ 0 ]->SetUpPointsOnLine(1, _2dMTLattice);
}

void
Lattice2d_mt :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = { P_f };
}

void
Lattice2d_mt :: initializeFrom(InputRecord &ir, int priority)
{
    ParameterManager &ppm =  this->giveDomain()->elementPPM;

    // first call parent
    LatticeTransportElement :: initializeFrom(ir, priority);
    PM_UPDATE_PARAMETER(dimension, ppm, ir, this->number, IPK_Lattice2d_mt_dim, priority) ;
    PM_UPDATE_PARAMETER(thickness, ppm, ir, this->number, IPK_Lattice2d_mt_thickness, priority) ;
    PM_UPDATE_PARAMETER(width, ppm, ir, this->number, IPK_Lattice2d_mt_width, priority) ;
    PM_UPDATE_PARAMETER(gpCoords, ppm, ir, this->number, IPK_Lattice2d_mt_gpcoords, priority) ;
    PM_UPDATE_PARAMETER(crackWidths.at(1), ppm, ir, this->number, IPK_Lattice2d_mt_crackwidth, priority) ;
    PM_UPDATE_PARAMETER(couplingFlag, ppm, ir, this->number, IPK_Lattice2d_mt_couplingflag, priority) ;
    PM_UPDATE_PARAMETER(couplingNumbers.at(1), ppm, ir, this->number, IPK_Lattice2d_mt_couplingnumber, priority) ;
}

void
Lattice2d_mt :: postInitialize()
{
    ParameterManager &ppm =  this->giveDomain()->elementPPM;
    this->LatticeTransportElement :: postInitialize();

    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_Lattice2d_mt_thickness) ;
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_Lattice2d_mt_width) ;
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_Lattice2d_mt_gpcoords) ;
    crackLengths.resize(1);
    crackLengths.at(1) = width;
    numberOfGaussPoints = 1;
}

double
Lattice2d_mt :: computeVolumeAround(GaussPoint *gp)
{
    return this->width * this->thickness * this->giveLength();
}

void
Lattice2d_mt :: computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    std :: unique_ptr< IntegrationRule > &iRule = integrationRulesArray [ 0 ];
    GaussPoint *gp = iRule->getIntegrationPoint(0);

    answer.resize(2, 2);
    answer.zero();
    answer.at(1, 1) = 1.;
    answer.at(1, 2) = -1.;
    answer.at(2, 1) = -1;
    answer.at(2, 2) = 1.;

    double k = static_cast< TransportMaterial * >( this->giveMaterial() )->giveCharacteristicValue(Conductivity, gp, tStep);
    double dV = this->computeVolumeAround(gp);
    double length = giveLength();
    double temp = k * dV / pow(length, 2.);
    answer.times(temp);
}


void
Lattice2d_mt :: computeCapacityMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    answer.resize(2, 2);
    answer.zero();
    answer.at(1, 1) = 2.;
    answer.at(1, 2) = 1.;
    answer.at(2, 1) = 1.;
    answer.at(2, 2) = 2.;
    double c = static_cast< TransportMaterial * >( this->giveMaterial() )->giveCharacteristicValue(Capacity, gp, tStep);
    double dV = this->computeVolumeAround(gp) / ( 6.0 * this->dimension );
    answer.times(c * dV);
}



void
Lattice2d_mt :: computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    std :: unique_ptr< IntegrationRule > &iRule = integrationRulesArray [ 0 ];

    FloatArray deltaX(3);
    FloatArray val, helpLoadVector, globalIPcoords;
    FloatMatrix nm;
    answer.clear();

    FloatArray gravityHelp(2);

    int nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        int n = bodyLoadArray.at(i);
        auto load  = domain->giveLoad(n);
        auto ltype = load->giveBCGeoType();

        if ( ltype == GravityPressureBGT ) {
            //Compute change of coordinates
            Node *nodeA = this->giveNode(1);
            Node *nodeB = this->giveNode(2);
            deltaX.at(1) = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
            deltaX.at(2) = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
            deltaX.at(3) = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);

            //Compute the local coordinate system
            GaussPoint *gp = iRule->getIntegrationPoint(0);

            gravityHelp.at(1) = 1.;
            gravityHelp.at(2) = -1.;

            double dV = this->computeVolumeAround(gp);
            load->computeValueAt(val, tStep, deltaX, mode);

            double k = static_cast< TransportMaterial * >( this->giveMaterial() )->giveCharacteristicValue(Conductivity, gp, tStep);

            double helpFactor = val.at(1) * k * dV;

            helpFactor /= pow(this->giveLength(), 2.);
            gravityHelp.times(helpFactor);

            if ( helpLoadVector.isEmpty() ) {
                helpLoadVector.resize(gravityHelp.giveSize() );
            }

            for ( int j = 1; j <= gravityHelp.giveSize(); j++ ) {
                helpLoadVector.at(j) += gravityHelp.at(j);
            }
        }

        answer.add(helpLoadVector);
    }
}


void
Lattice2d_mt :: giveGpCoordinates(FloatArray &answer)
{
    answer.resize(3);
    answer.at(1) = this->gpCoords.at(1);
    answer.at(2) = this->gpCoords.at(2);
}

int
Lattice2d_mt :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    answer.resize(3);
    answer.at(1) = this->gpCoords.at(1);
    answer.at(2) = this->gpCoords.at(2);

    return 1;
}

#define POINT_TOL 1.e-3

bool
Lattice2d_mt :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
{
    answer.resize(1);
    answer.at(1) = 0.;

    return true;
}


#ifdef __OOFEG


void
Lattice2d_mt :: drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
{
    OGC_PlotModeType mode = gc.giveIntVarPlotMode();

    if ( mode == OGC_rawGeometry ) {
        this->drawRawGeometry(gc, tStep);
        this->drawRawCrossSections(gc, tStep);
    } else {
        OOFEM_ERROR("unsupported mode");
    }
}




void Lattice2d_mt :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

void Lattice2d_mt :: drawRawCrossSections(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getCrossSectionColor() );
    EASValsSetLayer(OOFEG_RAW_CROSSSECTION_LAYER);

    FloatArray coords;
    this->giveCrossSectionCoordinates(coords);

    p [ 0 ].x = ( FPNum ) coords.at(1);
    p [ 0 ].y = ( FPNum ) coords.at(2);
    p [ 0 ].z = ( FPNum ) coords.at(3);
    p [ 1 ].x = ( FPNum ) coords.at(4);
    p [ 1 ].y = ( FPNum ) coords.at(5);
    p [ 1 ].z = ( FPNum ) coords.at(6);

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

void
Lattice2d_mt :: giveCrossSectionCoordinates(FloatArray &coords)
{
    double x1 = this->giveNode(1)->giveCoordinate(1);
    double y1 = this->giveNode(1)->giveCoordinate(2);
    double x2 = this->giveNode(2)->giveCoordinate(1);
    double y2 = this->giveNode(2)->giveCoordinate(2);

    //Compute normal and shear direction
    FloatArray normalDirection(2);
    FloatArray shearDirection(2);

    normalDirection.at(1) = x2 - x1;
    normalDirection.at(2) = y2 - y1;
    normalDirection.normalize();
    if ( normalDirection.at(2) == 0. ) {
        shearDirection.at(1) = 0.;
        shearDirection.at(2) = 1.;
    } else {
        shearDirection.at(1) = 1.0;
        shearDirection.at(2) = -normalDirection.at(1) / normalDirection.at(2);
    }

    shearDirection.normalize();

    coords.resize(6);
    coords.at(1) = this->gpCoords.at(1) - shearDirection.at(1) * this->width / 2.;
    coords.at(2) = this->gpCoords.at(2) - shearDirection.at(2) * this->width / 2.;
    coords.at(3) = 0.;
    coords.at(4) = this->gpCoords.at(1) + shearDirection.at(1) * this->width / 2.;
    coords.at(5) = this->gpCoords.at(2) + shearDirection.at(2) * this->width / 2.;
    coords.at(6) = 0.;
}

#endif
} // end namespace oofem
