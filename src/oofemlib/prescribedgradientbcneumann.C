/*
 * prescribedgradientbcneumann.C
 *
 *  Created on: Mar 5, 2014
 *      Author: svennine
 */

#include "prescribedgradientbcneumann.h"
#include "classfactory.h"
#include "node.h"
#include "masterdof.h"
#include "element.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "gausspoint.h"
#include "sparsemtrx.h"
#include "xfem/xfemelementinterface.h"
#include "xfem/integrationrules/discsegintegrationrule.h"

#include <cmath>

namespace oofem {
REGISTER_BoundaryCondition(PrescribedGradientBCNeumann);

PrescribedGradientBCNeumann::PrescribedGradientBCNeumann(int n, Domain * d):
PrescribedGradientBC(n, d)
{
    int nsd = d->giveNumberOfSpatialDimensions();
    int numComponents = 0;

    switch(nsd){
    case 1:
    	numComponents = 1;
    	break;

    case 2:
    	numComponents = 4;
    	break;

    case 3:
    	numComponents = 9;
    	break;
    }

    mSigmaIds.clear();
    mpSigmaHom = new Node(0, d); // Node number lacks meaning here.
    for ( int i = 0; i < numComponents; i++ ) {
        // Just putting in X_i id-items since they don't matter.
        int dofId = d->giveNextFreeDofID();
        mSigmaIds.followedBy(dofId);
    	mpSigmaHom->appendDof( new MasterDof( i + 1, mpSigmaHom, ( DofIDItem ) ( dofId ) ) );
    }

}

PrescribedGradientBCNeumann::~PrescribedGradientBCNeumann()
{
	if(mpSigmaHom != NULL) {
		delete mpSigmaHom;
		mpSigmaHom = NULL;
	}
}

DofManager* PrescribedGradientBCNeumann::giveInternalDofManager(int i)
{
	return mpSigmaHom;
}

void PrescribedGradientBCNeumann::scale(double s)
{

}

void PrescribedGradientBCNeumann::assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                            CharType type, ValueModeType mode,
                            const UnknownNumberingScheme &s, FloatArray *eNorm)
{
    Set *setPointer = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = setPointer->giveBoundaryList();

    IntArray loc, sigma_loc;  // For the displacements and stress respectively
    IntArray masterDofIDs, sigmaMasterDofIDs;
    mpSigmaHom->giveCompleteLocationArray(sigma_loc, s);

    if ( type == ExternalForcesVector ) {
        // The external forces have two contributions. On the additional equations for sigma, the load is simply the prescribed gradient.
        double rve_size = this->domainSize();
        FloatArray stressLoad;
        FloatArray gradVoigt;
        giveGradientVoigt(gradVoigt);

        stressLoad.beScaled(-rve_size, gradVoigt);
        answer.assemble(stressLoad, sigma_loc);

    } else if ( type == InternalForcesVector ) {
        FloatMatrix Ke;
        FloatArray fe_v, fe_s;
        FloatArray sigmaHom, e_u;

        // Fetch the current values of the stress;
        mpSigmaHom->giveCompleteUnknownVector(sigmaHom, mode, tStep);
        // and the master dof ids for sigmadev used for the internal norms
        mpSigmaHom->giveCompleteMasterDofIDArray(sigmaMasterDofIDs);

        // Assemble
        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            // Fetch the element information;
            e->giveLocationArray(loc, eid, s, &masterDofIDs);
            // Here, we could use only the nodes actually located on the boundary, but we don't.
            // Instead, we use all nodes belonging to the element, which is allowed because the
            // basis functions related to the interior nodes will be zero on the boundary.
            // Obviously, this is less efficient, so why do we want to do it this way?
            // Because it is easier when XFEM enrichments are present. /ES
            e->computeVectorOf(eid, mode, tStep, e_u);
            this->integrateTangent(Ke, e, boundary);

            // We just use the tangent, less duplicated code (the addition of sigmaDev is linear).
            fe_v.beProductOf(Ke, e_u);
            fe_s.beTProductOf(Ke, sigmaHom);

            // Note: The terms appear negative in the equations:
            fe_v.negated();
            fe_s.negated();

            answer.assemble(fe_s, loc); // Contributions to delta_v equations
            answer.assemble(fe_v, sigma_loc); // Contribution to delta_s_i equations
            if ( eNorm != NULL ) {
                eNorm->assembleSquared(fe_s, masterDofIDs);
                eNorm->assembleSquared(fe_v, sigmaMasterDofIDs);
            }
        }
    }
}

void PrescribedGradientBCNeumann::assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                      CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( eid == EID_MomentumBalance_ConservationEquation ) {
        eid = EID_MomentumBalance;
    }

    if ( eid != EID_MomentumBalance ) {
    	printf("returning.");
        return;
    }

    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == StiffnessMatrix || type == ElasticStiffnessMatrix ) {
        FloatMatrix Ke, KeT;
        IntArray loc_r, loc_c, sigma_loc_r, sigma_loc_c;
        Set *set = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = set->giveBoundaryList();

        // Fetch the columns/rows for the stress contributions;
        mpSigmaHom->giveCompleteLocationArray(sigma_loc_r, r_s);
        mpSigmaHom->giveCompleteLocationArray(sigma_loc_c, c_s);

        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            // Here, we could use only the nodes actually located on the boundary, but we don't.
            // Instead, we use all nodes belonging to the element, which is allowed because the
            // basis functions related to the interior nodes will be zero on the boundary.
            // Obviously, this is less efficient, so why do we want to do it this way?
            // Because it is easier when XFEM enrichments are present. /ES
            e->giveLocationArray(loc_r, eid, r_s);
            e->giveLocationArray(loc_c, eid, c_s);

            this->integrateTangent(Ke, e, boundary);
            Ke.negated();
            KeT.beTranspositionOf(Ke);

            answer->assemble(sigma_loc_r, loc_c, Ke); // Contribution to delta_s_i equations
            answer->assemble(loc_r, sigma_loc_c, KeT); // Contributions to delta_v equations
        }
    }
    else {
    	printf("Skipping assembly in PrescribedGradientBCNeumann::assemble().\n");
    }

}

void PrescribedGradientBCNeumann::giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, EquationID eid, CharType type,
                                const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    IntArray bNodes, dofids;
    IntArray loc_r, loc_c, sigma_loc_r, sigma_loc_c;
    int nsd = this->domain->giveNumberOfSpatialDimensions();
    DofIDItem id0 = this->domain->giveDofManager(1)->hasDofID(V_u) ? V_u : D_u; // Just check the first node if it has V_u or D_u.
    dofids.resize(nsd);
    for ( int i = 0; i < nsd; ++i ) {
        dofids(i) = id0 + i;
    }

    // Fetch the columns/rows for the stress contributions;
    mpSigmaHom->giveCompleteLocationArray(sigma_loc_r, r_s);
    mpSigmaHom->giveCompleteLocationArray(sigma_loc_c, c_s);

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    rows.resize( boundaries.giveSize() );
    cols.resize( boundaries.giveSize() );
    int i = 0;
    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );

        // Here, we could use only the nodes actually located on the boundary, but we don't.
        // Instead, we use all nodes belonging to the element, which is allowed because the
        // basis functions related to the interior nodes will be zero on the boundary.
        // Obviously, this is less efficient, so why do we want to do it this way?
        // Because it is easier when XFEM enrichments are present. /ES
        e->giveLocationArray(loc_r, eid, r_s);
        e->giveLocationArray(loc_c, eid, c_s);

        // For most uses, loc_r == loc_c, and sigma_loc_r == sigma_loc_c.
        rows [ i ] = loc_r;
        cols [ i ] = sigma_loc_c;
        i++;
        // and the symmetric part (usually the transpose of above)
        rows [ i ] = sigma_loc_r;
        cols [ i ] = loc_c;
        i++;
    }
}

void PrescribedGradientBCNeumann::computeField(FloatArray &sigma, EquationID eid, TimeStep *tStep)
{
    mpSigmaHom->giveUnknownVector(sigma, mSigmaIds, VM_Total, tStep);
}

void PrescribedGradientBCNeumann::integrateTangent(FloatMatrix &oTangent, Element *e, int iBndIndex)
{
    FloatArray normal, n;
    FloatMatrix nMatrix, E_n;
    FloatMatrix contrib;

    Domain *domain = e->giveDomain();
    XfemElementInterface *xfemElInt = dynamic_cast<XfemElementInterface*> (e);

    FEInterpolation *interp = e->giveInterpolation(); // Geometry interpolation

    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();

    // Interpolation order
    int order = interp->giveInterpolationOrder();
    IntegrationRule *ir = NULL;

    IntArray edgeNodes;
    FEInterpolation2d *interp2d = dynamic_cast<FEInterpolation2d*> (interp);
    if(interp2d == NULL) {
    	OOFEM_ERROR("failed to cast to FEInterpolation2d.")
    }
    interp2d->computeLocalEdgeMapping(edgeNodes, iBndIndex);

    const FloatArray &xS = *(e->giveDofManager(edgeNodes.at(1))->giveCoordinates());
    const FloatArray &xE = *(e->giveDofManager(edgeNodes.at(edgeNodes.giveSize()))->giveCoordinates());

    if(xfemElInt != NULL && domain->hasXfemManager() ) {
    	std::vector<Line> segments;
    	std::vector<FloatArray> intersecPoints;
    	xfemElInt->partitionEdgeSegment( iBndIndex, segments, intersecPoints );
        MaterialMode matMode = e->giveMaterialMode();
    	ir = new DiscontinuousSegmentIntegrationRule(1, e, segments, xS, xE);
    	int numPointsPerSeg = 1;
    	ir->SetUpPointsOnLine(numPointsPerSeg, matMode);
    }
    else {
    	ir = interp->giveBoundaryIntegrationRule(order, iBndIndex);
    }

    oTangent.clear();

    for ( GaussPoint *gp: *ir ) {
        FloatArray &lcoords = * gp->giveLocalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        // Evaluate the normal;
        double detJ = interp->boundaryEvalNormal(normal, iBndIndex, lcoords, cellgeo);

        interp->boundaryEvalN(n, iBndIndex, lcoords, cellgeo);
        // If cracks cross the edge, special treatment is necessary.
        // Exploit the XfemElementInterface to minimize duplication of code.
        if(xfemElInt != NULL && domain->hasXfemManager()) {
            // Compute global coordinates of Gauss point
            FloatArray globalCoord;

            interp->boundaryLocal2Global(globalCoord, iBndIndex, lcoords, cellgeo);

            // Compute local coordinates on the element
            FloatArray locCoord;
            e->computeLocalCoordinates(locCoord, globalCoord);

            xfemElInt->XfemElementInterface_createEnrNmatrixAt(nMatrix, locCoord, *e);
        } else {
            // Evaluate the velocity/displacement coefficients
            nMatrix.beNMatrixOf(n, nsd);
        }

        if ( nsd == 3 ) {
        	OOFEM_ERROR("not implemented for nsd == 3.")
        } else if ( nsd == 2 ) {

            E_n.resize(4, 2);
            E_n.at(1, 1) = normal.at(1);
            E_n.at(1, 2) = 0.0;

            E_n.at(2, 1) = 0.0;
            E_n.at(2, 2) = normal.at(2);

            E_n.at(3, 1) = normal.at(2);
            E_n.at(3, 2) = 0.0;

            E_n.at(4, 1) = 0.0;
            E_n.at(4, 2) = normal.at(1);


        } else {
            E_n.clear();
        }

        contrib.beProductOf(E_n, nMatrix);

        oTangent.add(detJ * gp->giveWeight(), contrib);
    }
    delete ir;
}


} /* namespace oofem */
