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

#include "transversereinfconstraint.h"
#include "classfactory.h"
#include "node.h"
#include "masterdof.h"
#include "element.h"
#include "feinterpol.h"
#include "gausspoint.h"
#include "sparsemtrx.h"
#include "xfem/xfemelementinterface.h"
#include "xfem/integrationrules/discsegintegrationrule.h"
#include "timestep.h"
#include "function.h"
#include "sparselinsystemnm.h"
#include "unknownnumberingscheme.h"
#include "engngm.h"
#include "mathfem.h"

namespace oofem {
REGISTER_BoundaryCondition(TransverseReinfConstraint);

TransverseReinfConstraint :: TransverseReinfConstraint(int n, Domain *d) :
    ActiveBoundaryCondition(n, d),
    lmLambda( new Node(0, d) )
{
    int dofId = d->giveNextFreeDofID();
    lmLambdaIds.followedBy(dofId);
    lmLambda->appendDof( new MasterDof( lmLambda.get(), ( DofIDItem ) ( dofId ) ) );
}

TransverseReinfConstraint :: ~TransverseReinfConstraint()
{
}


void TransverseReinfConstraint :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, steelElSet, _IFT_TransverseReinfConstraint_SteelElSet);
    IR_GIVE_FIELD(ir, conElBoundSet, _IFT_TransverseReinfConstraint_ConElBoundSet);

    return ActiveBoundaryCondition :: initializeFrom(ir);
}


void TransverseReinfConstraint :: giveInputRecord(DynamicInputRecord &input)
{
    ActiveBoundaryCondition :: giveInputRecord(input);
}


DofManager *TransverseReinfConstraint :: giveInternalDofManager(int i)
{
    return lmLambda.get();
}


void TransverseReinfConstraint :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                 CharType type, ValueModeType mode,
                                                 const UnknownNumberingScheme &s,
                                                 FloatArray *eNorm,
                                                 void* lock)
{
    IntArray lambda_loc;
    lmLambda->giveLocationArray(lmLambdaIds, lambda_loc, s);

    if ( type == ExternalForcesVector ) {
        // The external forces have two contributions. On the additional equations for lambda, the load is simply zero.
        FloatArray stressLoad = Vec1(0);
        answer.assemble(stressLoad, lambda_loc);
    } else if ( type == InternalForcesVector ) {
        FloatMatrix Klam;
        FloatArray fe_lam, fe_u;
        FloatArray lambda, u_s, u_c, u;
        IntArray loc, loc_s, loc_c, masterDofIDs;

        // Fetch the current values of the Lagrange multiplier;
        lmLambda->giveUnknownVector(lambda, lmLambdaIds, mode, tStep);

        // Assemble
        Set *steelSet = this->giveDomain()->giveSet(steelElSet);
        const IntArray &elements = steelSet->giveElementList();

        Set *concreteSet = this->giveDomain()->giveSet(conElBoundSet);
        const IntArray &boundaries = concreteSet->giveBoundaryList();

        if ( boundaries.giveSize() == elements.giveSize()*2 ) {
            for ( int pos = 1; pos <= elements.giveSize(); ++pos) {
                Element *es = this->giveDomain()->giveElement( elements.at(pos) );
                Element *ec = this->giveDomain()->giveElement( boundaries.at(2*pos - 1) );
                int boundary = boundaries.at(2*pos);

                // Fetch the element information;
                es->giveLocationArray(loc_s, s);

                // Fetch the nodal displacements
                // since computeVectorOf gives the unknown vector in local coordinate system, it needs to be rotated to global coordinates
                FloatMatrix G2L;
                es->computeGtoLRotationMatrix(G2L);
                es->computeVectorOf(mode, tStep, u_s);
                u_s.rotatedWith(G2L, 't');

                //Fetch boundary information
                IntArray bndNodes;
                bndNodes = ec->giveBoundaryEdgeNodes(boundary);
                ec->giveBoundaryLocationArray(loc_c, bndNodes, s, &masterDofIDs);
                masterDofIDs.resizeWithValues(2); //hack to get corresponding masterDofIDs in computeBoundaryVectorOf
                ec->computeBoundaryVectorOf(bndNodes, masterDofIDs, mode, tStep, u_c);

                //Assemble the location and u arrays
                loc=loc_s;
                loc.followedBy(loc_c);

                u=FloatArray::fromConcatenated({u_s,u_c});

                //Compute the stiffness matrix expansion
                this->integrateTangent(Klam, es, ec, boundary);

                //Compute the contribution to internal force vector
                fe_lam.beProductOf(Klam, u);
                fe_u.beTProductOf(Klam, lambda);

                //Negate the terms for symmetry
                fe_lam.negated();
                fe_u.negated();

                //Assemble
                answer.assemble(fe_u, loc);
                answer.assemble(fe_lam, lambda_loc);
            }
        } else {
            OOFEM_ERROR("Steel/concrete element mismatch")
        }
    }
}

void TransverseReinfConstraint :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                           CharType type,
                                           const UnknownNumberingScheme &r_s,
                                           const UnknownNumberingScheme &c_s,
                                           double scale,
                                           void* lock)
{
    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == ElasticStiffnessMatrix ) {
        FloatMatrix Klam, KlamT;
        IntArray loc_r, loc_c, lambda_loc_r, lambda_loc_c;
        IntArray loc_us_r, loc_us_c, loc_uc_r, loc_uc_c;
        IntArray masterDofIDs;

        //Fetch the rows/columns for the Lagrange multiplier contributions
        lmLambda->giveLocationArray(lmLambdaIds, lambda_loc_r, r_s);
        lmLambda->giveLocationArray(lmLambdaIds, lambda_loc_c, c_s);

        //Contributions from steel
        Set *steelSet = this->giveDomain()->giveSet(steelElSet);
        const IntArray &elements = steelSet->giveElementList();

        //Contributions from concrete
        Set *concreteSet = this->giveDomain()->giveSet(conElBoundSet);
        const IntArray &boundaries = concreteSet->giveBoundaryList();

        if ( boundaries.giveSize() == elements.giveSize()*2 ) {
            for ( int pos = 1; pos <= elements.giveSize(); ++pos) {
                Element *es = this->giveDomain()->giveElement( elements.at(pos) );
                Element *ec = this->giveDomain()->giveElement( boundaries.at(2*pos - 1) );
                int boundary = boundaries.at(2*pos);

                // Fetch the element information;
                es->giveLocationArray(loc_us_r, r_s);
                es->giveLocationArray(loc_us_c, c_s);

                //Fetch boundary information
                IntArray bndNodes;
                bndNodes = ec->giveBoundaryEdgeNodes(boundary);
                ec->giveBoundaryLocationArray(loc_uc_r, bndNodes, r_s, &masterDofIDs);
                ec->giveBoundaryLocationArray(loc_uc_c, bndNodes, c_s, &masterDofIDs);

                //Assemble the location arrays
                loc_c=loc_us_c;
                loc_c.followedBy(loc_uc_c);
                loc_r=loc_us_r;
                loc_r.followedBy(loc_uc_r);

                //Compute the stiffness matrix expansion
                this->integrateTangent(Klam, es, ec, boundary);

                Klam.negated();
                KlamT.beTranspositionOf(Klam);

                answer.assemble(lambda_loc_r, loc_c, Klam);
                answer.assemble(loc_r, lambda_loc_c, KlamT);
            }
        } else {
            OOFEM_ERROR("Steel/concrete element mismatch")
        }
    } else {
        OOFEM_LOG_DEBUG("Skipping assembly in TransverseReinfConstraint::assemble().\n");
    }
}

void TransverseReinfConstraint :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                                       const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    IntArray loc_us_r, loc_us_c, loc_uc_r, loc_uc_c;
    IntArray loc_r, loc_c, lambda_loc_r, lambda_loc_c;
    IntArray masterDofIDs;

    // Fetch the columns/rows for the stress contributions;
    lmLambda->giveLocationArray(lmLambdaIds, lambda_loc_r, r_s);
    lmLambda->giveLocationArray(lmLambdaIds, lambda_loc_c, c_s);

    //Contributions from steel
    Set *steelSet = this->giveDomain()->giveSet(steelElSet);
    const IntArray &selements = steelSet->giveElementList();

    //Contributions from concrete
    Set *concreteSet = this->giveDomain()->giveSet(conElBoundSet);
    const IntArray &boundaries = concreteSet->giveBoundaryList();

    if ( boundaries.giveSize() == selements.giveSize()*2 ) {
        rows.resize( boundaries.giveSize() + 2*selements.giveSize() );
        cols.resize( boundaries.giveSize() + 2*selements.giveSize() );

        int i=0;
        for ( int pos = 1; pos <= selements.giveSize(); ++pos) {
            Element *es = this->giveDomain()->giveElement( selements.at(pos) );
            Element *ec = this->giveDomain()->giveElement( boundaries.at(2*pos - 1) );
            int boundary = boundaries.at(2*pos);

            // Fetch the element information;
            es->giveLocationArray(loc_us_r, r_s);
            es->giveLocationArray(loc_us_c, c_s);

            //Fetch boundary information
            IntArray bndNodes;
            bndNodes = ec->giveBoundaryEdgeNodes(boundary);
            ec->giveBoundaryLocationArray(loc_uc_r, bndNodes, r_s, &masterDofIDs);
            ec->giveBoundaryLocationArray(loc_uc_c, bndNodes, c_s, &masterDofIDs);

            loc_c=loc_us_c;
            loc_c.followedBy(loc_uc_c);
            loc_r=loc_us_r;
            loc_r.followedBy(loc_uc_r);

            // For most uses, loc_us_r == loc_us_c, and lambda_loc_r == lambda_loc_c.
            rows [ i ] = loc_r;
            cols [ i ] = lambda_loc_c;
            i++;
            // and the symmetric part (usually the transpose of above)
            rows [ i ] = lambda_loc_r;
            cols [ i ] = loc_c;
            i++;
        }
    } else {
        OOFEM_ERROR("Steel/concrete element mismatch")
    }
}

void TransverseReinfConstraint :: computeField(FloatArray &lambda, TimeStep *tStep)
{
    lmLambda->giveUnknownVector(lambda, lmLambdaIds, VM_Total, tStep);
}



void TransverseReinfConstraint :: integrateTangent(FloatMatrix& oTangent, Element* es, Element* ec, int iBndIndex)
{
    FloatMatrix Kslam, Kclam;
    int nslam, nclam, nsd=1;

    //Compute contributions from steel and concrete
    this->integrateTangentOnSteel(Kslam, es);
    this->integrateTangentOnConcrete(Kclam, ec, iBndIndex);
    nslam = Kslam.giveNumberOfColumns();
    nclam = Kclam.giveNumberOfColumns();

    //Assemble the tangent contribution
    oTangent.resize(nsd,nslam + nclam);
    for (int i=1; i <= nsd; ++i) {
        for (int j=1; j <= nslam; ++j) {
            oTangent.at(i,j) = Kslam.at(i,j);
        }
        for (int j=nslam+1; j <= nslam+nclam; ++j) {
            oTangent.at(i,j) = Kclam.at(i,j-nslam);
        }
    }
}

void TransverseReinfConstraint :: integrateTangentOnSteel(FloatMatrix& oTangent, Element* e)
{
    //reminder: Kslam must be 1x6
    FloatArray normal;
    FloatMatrix Ns, contrib;
    IntArray DofMask;

    e->giveDofManDofIDMask(1, DofMask);
    int ndof = DofMask.giveSize();

    FEInterpolation *interp = e->giveInterpolation();
    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir;
    ir = interp->giveIntegrationRule(order, e->giveGeometryType());

    oTangent.clear();

    for (auto &gp : *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        FloatArray nu;
        interp->evalN(nu, lcoords, cellgeo);
        double detJ = interp->boundaryEvalNormal(normal, 1, lcoords, cellgeo);

        //Constant traction interpolation
        FloatMatrix Nlam(ndof,1);
        // Note: this works only if the reinforcement is horizontal/vertical!
        // Regardless of interpolation and direction of normal vector taking an absolute value gives
        // consistent direction of Lagrange multiplier in this case.
        // Otherwise the order of nodes for the line elements would make a difference (inward/outward pointing normal).
        // Ideally, it shouldn't matter in which direction (bottom->up, top->down) the reinforcement elements are defined.
        Nlam(0,0) = fabs(normal.at(1)); //x component of LM
        Nlam(1,0) = fabs(normal.at(2)); //y component of LM
        Nlam(2,0) = 0.; //rotation component of LM

        Ns.beNMatrixOf(nu, ndof);

        contrib.beTProductOf(Nlam,Ns);
        oTangent.add(detJ * gp->giveWeight(), contrib);
    }
}

void TransverseReinfConstraint :: integrateTangentOnConcrete(FloatMatrix &oTangent, Element *e, int iBndIndex)
{
    //reminder: Kclam must be 1x4
    FloatArray normal;
    FloatMatrix Nc, contrib;
    IntArray DofMask;

    e->giveDofManDofIDMask(1, DofMask);
    int ndof = DofMask.giveSize();

    FEInterpolation *interp = e->giveInterpolation();
    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir;
    ir = interp->giveBoundaryIntegrationRule(order, iBndIndex, e->giveGeometryType());

    oTangent.clear();

    for (auto &gp : *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        FloatArray nu;
        interp->boundaryEvalN(nu, iBndIndex, lcoords, cellgeo);
        double detJ = interp->boundaryEvalNormal(normal, iBndIndex, lcoords, cellgeo);

        //Constant traction interpolation
        FloatMatrix Nlam(ndof,1);
        // Note: this works only if the reinforcement is horizontal/vertical!
        // Regardless of interpolation and direction of normal vector taking an absolute value gives
        // consistent direction of Lagrange multiplier in this case.
        // Otherwise it would make a difference on which side of the reinforcement bar the solid elements are located,
        // but it shouldn't matter.
        Nlam(0,0) = fabs(normal.at(1)); //x component of LM
        Nlam(1,0) = fabs(normal.at(2)); //y component of LM

        Nc.beNMatrixOf(nu, ndof);

        contrib.beTProductOf(Nlam,Nc);
        oTangent.add(detJ * gp->giveWeight(), contrib);
        oTangent.negated(); //this term is negative in the formulation
    }
}
} /* namespace oofem */
