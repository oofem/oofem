#include "node2segmentlagrangemultipliercontact.h"

namespace oofem {
    REGISTER_BoundaryCondition(Node2SegmentLagrangianMultiplierContact);



    IRResultType Node2SegmentLagrangianMultiplierContact::initializeFrom(InputRecord * ir)
    {
        IRResultType result;
        this->useTangent = ir->hasField(_IFT_Node2SegmentLagrangianMultiplierContact_useTangent);

        IR_GIVE_FIELD(ir, this->segmentSet, _IFT_Node2SegmentLagrangianMultiplierContact_segmentSet);
        IR_GIVE_FIELD(ir, this->nodeSet, _IFT_Node2SegmentLagrangianMultiplierContact_nodeSet);

        this->lm_num = segmentSet.giveSize() * nodeSet.giveSize();

        // these are internal lagrange multipliers used to enforce the no penetration for the contacting bodies
        for ( int pos = 0; pos < lm_num; ++pos ) {
            lmdm.emplace_back(/*std::make_unique<Node>*/ new Node(0, domain));
            lmdm.at(pos)->appendDof(new MasterDof(this->lmdm.at(pos), (DofIDItem)(this->giveDomain()->giveNextFreeDofID())));
        }

        IR_GIVE_OPTIONAL_FIELD(ir, this->prescribedNormal, _IFT_Node2SegmentLagrangianMultiplierContact_prescribedNormal);


        return ActiveBoundaryCondition::initializeFrom(ir);
    }

    void Node2SegmentLagrangianMultiplierContact::assemble(SparseMtrx & answer, TimeStep * tStep, CharType type, const UnknownNumberingScheme & r_s, const UnknownNumberingScheme & c_s, double scale)
    {
        if ( !this->useTangent || type != TangentStiffnessMatrix ) {
            return;
        }

        FloatMatrix K;
        IntArray loc, n_loc;

        IntArray dofIdArray = giveDomain()->giveDefaultNodeDofIDArry();
        /*IntArray dofIdArray = {
            D_u, D_v
        };*/

        std::vector< IntArray >lambdaeq;
        this->giveLagrangianMultiplierLocationArray(r_s, lambdaeq);
        int lmpos = 1;
        //iterate over all pairs of nodes and segments
        for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); ++nodePos ) {
            for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {

                Node* node = this->giveDomain()->giveNode(nodeSet.at(nodePos));
                ContactSegment* segment = (this->giveDomain()->giveContactSegment(segmentSet.at(segmentPos)));

                //node follows segment because that is how ContactSegments do it
                segment->giveLocationArray(dofIdArray, loc, r_s);
                node->giveLocationArray(dofIdArray, n_loc, c_s);
                loc.followedBy(n_loc);

                double gap = this->computeTangentFromContact(K, node, segment, tStep);
                if ( gap >= 0 ) {
                    // to make the equation system regular in the case there is no contact we initialize the allocated equation to the following form 1*labmda = 0, forcing lagrange multiplier of inactive condition to be zero.
                    FloatArray one(1);
                    one.at(1) = 1;
                    answer.assemble(lambdaeq.at(lmpos - 1), one);
                }
                else {
                    answer.assemble(loc, lambdaeq.at(lmpos - 1), K);
                    FloatMatrix Kt;
                    Kt.beTranspositionOf(K);
                    answer.assemble(lambdaeq.at(lmpos - 1), loc, Kt);
                }
                lmpos++;
            }
        }
    }

    void Node2SegmentLagrangianMultiplierContact::assembleVector(FloatArray & answer, TimeStep * tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme & s, FloatArray * eNorms)
    {
        IntArray dofIdArray = giveDomain()->giveDefaultNodeDofIDArry();
        /*IntArray dofIdArray = {
            D_u, D_v
        };*/

        if ( type == InternalForcesVector ) {
            // assemble lagrangian multiplier contribution to residuum
            // assemble location array
            std::vector< IntArray >lambdaeq;
            IntArray loc, n_loc;
            FloatArray n, fext;
            int lmpos = 1;

            this->giveLagrangianMultiplierLocationArray(s, lambdaeq);

            for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); ++nodePos ) {
                for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {

                    Node* node = this->giveDomain()->giveNode(nodeSet.at(nodePos));
                    ContactSegment* segment = (this->giveDomain()->giveContactSegment(segmentSet.at(segmentPos)));

                    this->computeNvMatrixAt(n, node, segment, tStep);
                    Dof *mdof = *(lmdm.at(lmpos - 1)->begin());

                    n.times(mdof->giveUnknown(mode, tStep));

                    this->computeExternalForcesFromContact(fext, node, segment, tStep);

                    segment->giveLocationArray(dofIdArray, loc, s);
                    node->giveLocationArray(dofIdArray, n_loc, s);
                    loc.followedBy(n_loc);

                    answer.assemble(n, loc);
                    answer.assemble(fext, lambdaeq.at(lmpos - 1));

                    lmpos++;
                }
            }
        }
    }

    void Node2SegmentLagrangianMultiplierContact::giveLocationArrays(std::vector<IntArray>& rows, std::vector<IntArray>& cols, CharType type, const UnknownNumberingScheme & r_s, const UnknownNumberingScheme & c_s)
    {
        IntArray r_loc, c_loc;
        rows.resize(3 * lm_num);
        cols.resize(3 * lm_num);
        IntArray dofIdArray = giveDomain()->giveDefaultNodeDofIDArry();
        /*IntArray dofIdArray = {
            D_u, D_v
        };*/

        std::vector< IntArray >lambdaeq;
        this->giveLagrangianMultiplierLocationArray(r_s, lambdaeq);

        int lmpos = 1;
        for ( int nodePos = 1; nodePos <= nodeSet.giveSize(); ++nodePos ) {
            for ( int segmentPos = 1; segmentPos <= segmentSet.giveSize(); segmentPos++ ) {
                Node* node = this->giveDomain()->giveNode(nodeSet.at(nodePos));
                ContactSegment* segment = (this->giveDomain()->giveContactSegment(segmentSet.at(segmentPos)));

                segment->giveLocationArrays(dofIdArray, r_loc, r_s);
                node->giveLocationArray(dofIdArray, c_loc, c_s);

                // column block
                rows[0 + 3 * (lmpos - 1)] = r_loc;
                cols[0 + 3 * (lmpos - 1)] = lambdaeq.at(lmpos - 1);
                // row block
                rows[1 + 3 * (lmpos - 1)] = c_loc;
                cols[1 + 3 * (lmpos - 1)] = lambdaeq.at(lmpos - 1);
                // diagonal enry (some sparse mtrx implementation requaire this)
                rows[2 + 3 * (lmpos - 1)] = lambdaeq.at(lmpos - 1);
                cols[2 + 3 * (lmpos - 1)] = lambdaeq.at(lmpos - 1);

                lmpos++;
            }
        }
    }

    void Node2SegmentLagrangianMultiplierContact::giveLagrangianMultiplierLocationArray(const UnknownNumberingScheme & r_s, std::vector<IntArray>& answer)
    {
        answer.resize(lm_num);
        IntArray dofIdArray = giveDomain()->giveDefaultNodeDofIDArry();
        /*IntArray dofIdArray = {
            D_u, D_v
        };*/

        // assemble location array
        IntArray l(1);
        for ( int i = 0; i < lm_num; i++ ) {
            l.at(1) = r_s.giveDofEquationNumber(*lmdm.at(i)->begin());
            answer.at(i) = l;
        }
    }

    double Node2SegmentLagrangianMultiplierContact::computeTangentFromContact(FloatMatrix & answer, Node * node, ContactSegment * segment, TimeStep * tStep)
    {
        double gap;
        FloatArray Nv;
        this->computeGap(gap, node, segment, tStep);
        this->computeNvMatrixAt(Nv, node, segment, tStep);
        answer.initFromVector(Nv, false);

        return gap;
    }

    void Node2SegmentLagrangianMultiplierContact::computeGap(double & answer, Node * node, ContactSegment * segment, TimeStep * tStep)
    {
        answer = segment->computePenetration(node, tStep);
    }

    void Node2SegmentLagrangianMultiplierContact::computeNvMatrixAt(FloatArray & answer, Node * node, ContactSegment * segment, TimeStep * tStep)
    {
        ////computeNormal is expected to return an integrated term
        //// int across seg (N^T), where N = element Nmatrix (extended by zeros for the node)

        //FloatArray normal;
        //FloatMatrix extendedN, extendedNTranspose;

        //if ( prescribedNormal.giveSize() == node->giveNumberOfDofs() ) {
        //    normal = prescribedNormal;
        //}
        //else {
        //    segment->computeNormal(normal, node, tStep);
        //    double norm = normal.computeNorm();
        //    if ( norm != 0 ) {
        //        normal.times(1. / norm);
        //    }
        //}

        //segment->computeSegmentNMatrix(extendedN, node, tStep);
        ////normal should be given just as N^t * n;
        //answer.beTProductOf(extendedN, normal);

		FloatArray normal;
		FloatMatrix segmentN, extendedN;
		int ndof = node->giveNumberOfDofs();

		if (prescribedNormal.giveSize() == node->giveNumberOfDofs()) {
			normal = prescribedNormal;
		}
		else {
			segment->computeNormal(normal, node, tStep);
		}
		int norm = normal.computeNorm();
		if (norm) {
			normal.normalize();
		}

		segment->computeSegmentNMatrix(segmentN, node, tStep);

		if (segmentN.giveNumberOfRows() != ndof) {
			OOFEM_ERROR("Dimension mismatch between node and contact segment");
		}

		//append extension. Currenty appended to the END (in penalty/large strain formulation, it is the other way around)
		extendedN.resize(segmentN.giveNumberOfRows(), segmentN.giveNumberOfColumns() + 2);
		extendedN.zero();
		FloatMatrix extension(ndof, ndof);
		extension.beUnitMatrix();
		extension.times(-1); //extension is negated, unlike penalty, where segmentN is negated instead
		extendedN.setSubMatrix(extension, 1, segmentN.giveNumberOfColumns() + 1);
		extendedN.setSubMatrix(segmentN, 1, 1);

		//Nv should be given just as N^t * n;
		answer.beTProductOf(extendedN, normal);
    }

    void Node2SegmentLagrangianMultiplierContact::computeExternalForcesFromContact(FloatArray & answer, Node * node, ContactSegment * segment, TimeStep * tStep)
    {
        answer.resize(1);
        this->computeGap(answer.at(1), node, segment, tStep);
        if ( answer.at(1) >= 0.0 ) {
            answer.at(1) = 0.0;
        }
        answer.times(-1.);
    }

}//end namespace oofem
