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


#include "sm/Contact/structuralpenaltycontactbc.h"
#include "domain.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "classfactory.h"
#include "mathfem.h"
#include <ranges>
#include <vector>
#include "floatarray.h"
#include "Contact/contactsearch.h"
#include "Contact/contactsearchsweepandprune.h"
#include "timestep.h"
namespace oofem {
REGISTER_BoundaryCondition(StructuralPenaltyContactBoundaryCondition);


/** temporary auxiliary function */
bool frictionShouldBeConsidered(double friction, TimeStep *tStep) {
	return friction > 0 && tStep->giveNumber() > 1;
}



void StructuralPenaltyContactBoundaryCondition::computeTangentFromContact(FloatMatrix &answer, ContactPair *contactPair, TimeStep *tStep)
{
    answer.zero();
    auto normal_gap = contactPair->giveNormalGap();
    if ( normal_gap <= 1.e-15 ) {
      contactPair->initContactPoint();
      auto normal = contactPair->giveNormalVector();
      auto area =  normal.computeNorm();
      normal /= area;
      // get the contact element shape function values at this node
      FloatMatrix N;
      contactPair->computeNmatrix(N);
      //
      double normal_traction;
      FloatArray tangential_traction, tangential_traction_trial;
      ContactProcess mode = ContactProcess::None;
      this->computeTractions(normal_traction, tangential_traction, tangential_traction_trial, mode, contactPair, tStep);
	  double abs_normal_traction = std::abs(normal_traction);
      /****************** Normal part of the stiffness matrix ******************************/
      //@todo: update
      double dA = 1;//area;
      //auto dA = contactPair->giveArea();
      // Equation(7.17) page 189
      answer.plus_Nt_a_otimes_b_B(N, normal, normal, N,  penalty_normal * dA);
      //
      auto tangent_vectors = contactPair->giveTangentVectors();
      auto contravariant_metric = this->computeContravariantMetric(tangent_vectors);
      // calculate curvature tensor kappa
      FloatMatrix curvature;
      contactPair->computeCurvature(curvature, tStep);  
      std::vector<FloatMatrix> dNs;
      contactPair->compute_dNdxi_matrices(dNs);
      FloatMatrix k_rot, k_curv;
      for(int i = 0; i < surface_dimension; i++) {
		auto tangent_i = tangent_vectors[i];
	    auto dNi = dNs[i];
		//for (auto && [i, tangent_i] : std::views::enumerate(tangent_vectors)) {
		for(int j = 0; j < surface_dimension; j++) {
		  auto dNj = dNs[j];
	  	  //for (auto && [j, dNj] : std::ranges::views::enumerate(Bs)) {
	  	  auto tangent_j = tangent_vectors[j];
		  double a_ij = contravariant_metric(i, j);
		  double h_ij = curvature(i, j);
	      // Equation(7.18)a page 189
	      k_rot.plus_Nt_a_otimes_b_B(dNj, normal, tangent_i, N,  a_ij * normal_traction * dA);
	      // Equation(7.18)b page 189
	      k_rot.plus_Nt_a_otimes_b_B(N, tangent_j, normal, dNi,   a_ij * normal_traction * dA);
	      // Equation(7.19) page 189
	      k_curv.plus_Nt_a_otimes_b_B(N, tangent_i, tangent_j, N,  h_ij * normal_traction * dA);
		}
      }
      answer +=  k_rot + k_curv;
      /****************** End of the Normal part of the stiffness matrix ******************************/
      if (frictionShouldBeConsidered(friction, tStep)) {
	    /****************** Frictional part of the stiffness matrix ******************************/
		tangential_traction = tangential_traction_trial;
	    if(mode == ContactProcess::Sticking) {
	      /****  STICK ****/
	      FloatMatrix k_st_m, k_st_r, k_st_c;
		  double factor;
	      /*for (auto && [i, rho_i] : std::ranges::views::enumerate(tangent_vectors)) {
	        for (auto && [j, rho_j] : std::ranges::views::enumerate(tangent_vectors)) {*/
	      for (int i = 0; i < this->surface_dimension; i++) {
	        auto rho_i = tangent_vectors[i];
			double t_i = tangential_traction(i);
	        for (int j = 0; j < this->surface_dimension; j++) {
	          auto rho_j = tangent_vectors[j];
			  double a_ij = contravariant_metric(i, j);
	          const auto& dNj = dNs[j];
			  double h_ij = curvature(i, j);
	          // Equation(7.26), page 192  
			  factor = - penalty_tangential * a_ij;
	          k_st_m.plus_Nt_a_otimes_b_B(N, rho_i, rho_j, N, factor);
	          // Equation (7.28), page 192
			  factor = + t_i * h_ij;
	          k_st_c.plus_Nt_a_otimes_b_B( N, rho_j, normal, N, factor);
	          k_st_c.plus_Nt_a_otimes_b_B( N, normal, rho_j, N, factor);
			  //for (auto && [k, tangent_k] : std::ranges::views::enumerate(tangent_vectors)) {
	          for (int k = 0; k < surface_dimension; k++) {
				auto rho_k = tangent_vectors[k];
				double a_ik = contravariant_metric(i, k);
				double a_jk = contravariant_metric(j, k);
				//for (auto && [l, tangent_l] : std::ranges::views::enumerate(tangent_vectors)) {
				for (int l = 0; l < surface_dimension; l++) {
				  auto rho_l = tangent_vectors[l];
				  double a_il = contravariant_metric(i, l);
				  double a_jl = contravariant_metric(j, l);
				  // Equation (7.27), page 192
				  factor = - t_i * a_il * a_jk;
				  k_st_r.plus_Nt_a_otimes_b_B( N, rho_k, rho_l, dNj, factor);
				  factor = - t_i * a_ik * a_jl;
				  k_st_r.plus_Nt_a_otimes_b_B( dNj, rho_k, rho_l, N, factor);
				}
	          }
	        }
	      }
		  FloatMatrix k_stick;
		  k_stick += k_st_m;
		  k_stick += k_st_r;
		  k_stick += k_st_c;
		  k_stick.times(dA);
	      answer += k_stick;
	      /****  END OF STICK ****/
	    } else {
	      /****  SLIP ****/
	      // see table 8.13 on the papge 252
	      // Equation (7.35), page 194
	      FloatMatrix k_sl_constitutive_non_symmetric; // (7.35a), page 194
	      FloatMatrix k_sl_constitutive_symmetric_1;   // (7.35b), page 194
	      FloatMatrix k_sl_constitutive_symmetric_2;   // (7.35c), page 194
	      FloatMatrix k_sl_rotational_symmetric;       // (7.35d), page 194
	      FloatMatrix k_sl_curvature_symmetric;        // (7.35e), page 194
	      FloatMatrix k_sl_curvature_non_symmetric;    // (7.35f), page 194
	      // Equation (7.35a), page 194
		  double t_norm_squared = 0.0;
		  for (int i = 0; i < surface_dimension; i++) {
			double t_i = tangential_traction(i);
			for (int j = 0; j < surface_dimension; j++) {
			  double t_j = tangential_traction(j);
			  double a_ij = contravariant_metric(i, j);
			  t_norm_squared += t_i * t_j * a_ij;
			}
		  }
	      const double t_norm = std::sqrt(t_norm_squared);
	      const double t_norm3 = std::pow(t_norm,3);
	      double factor;
	      for (int i = 0; i < surface_dimension; i++) {
			const auto& rho_i = tangent_vectors[i];
			const double t_i = tangential_traction(i);
			for (int j = 0; j < surface_dimension; j++) {
			  const auto& rho_j = tangent_vectors[j];
			  const double t_j = tangential_traction(j);
			  const double a_ij = contravariant_metric(i, j);
			  const auto& dNj = dNs[j];
			  const double h_ij = curvature(i,j);
			  // Equation (7.35a), page 194
			  factor = -penalty_normal * friction * t_i / t_norm * a_ij;
			  k_sl_constitutive_non_symmetric.plus_Nt_a_otimes_b_B(N, rho_j, normal, N, factor);
			  // Equation (7.35b), page 194
			  factor = -penalty_tangential * friction * abs_normal_traction / t_norm * a_ij;
			  k_sl_constitutive_symmetric_1.plus_Nt_a_otimes_b_B(N, rho_i, rho_j, N, factor);
			  // Equation (7.35e), page 194
			  factor = +friction * abs_normal_traction * t_i / t_norm * h_ij;
			  k_sl_curvature_symmetric.plus_Nt_a_otimes_b_B(N, rho_j, normal, N, factor);
			  k_sl_curvature_symmetric.plus_Nt_a_otimes_b_B(N, normal, rho_j, N, factor);
			  for (int k = 0; k < surface_dimension; k++) {
			    const auto& rho_k = tangent_vectors[k];
			    const double a_ik = contravariant_metric(i, k);
			    const double a_jk = contravariant_metric(j, k);
			    for (int l = 0; l < surface_dimension; l++) {
			      const auto& rho_l = tangent_vectors[l];
			      const double a_il = contravariant_metric(i, l);
			      const double a_jl = contravariant_metric(j, l);
			      // Equation (7.35c), page 194
			      factor = +penalty_tangential * friction * abs_normal_traction * t_i * t_j / t_norm3 * a_ik * a_jl;
			      k_sl_constitutive_symmetric_2.plus_Nt_a_otimes_b_B(N, rho_k, rho_l, N, factor);
			      // Equation (7.35d), page 194
			      const double factorCommon = -friction * abs_normal_traction * t_i / t_norm;
			      factor = factorCommon * a_il * a_jk;
			      k_sl_rotational_symmetric.plus_Nt_a_otimes_b_B(N, rho_k, rho_l, dNj, factor);
			      factor = factorCommon * a_ik * a_jl;
			      k_sl_rotational_symmetric.plus_Nt_a_otimes_b_B(dNj, rho_k, rho_l, N, factor);
			      for (int m = 0; m <= surface_dimension; m++) {
			        // Equation (7.35f), page 194
			        // @todo
			      }
			    }
	    	  }
	        }
	      }
		  FloatMatrix k_slip;
		  k_slip += k_sl_constitutive_non_symmetric;
		  k_slip += k_sl_constitutive_symmetric_1;
		  k_slip += k_sl_constitutive_symmetric_2;
		  k_slip += k_sl_rotational_symmetric;
		  k_slip += k_sl_curvature_symmetric;
		  k_slip += k_sl_curvature_non_symmetric;
		  k_slip.times(dA);
	      answer += k_slip;
	    }
      }
    }
}

  
FloatMatrix
StructuralPenaltyContactBoundaryCondition ::computeCovariantMetric(const std::vector<FloatArray> &tangent_vectors)
{
  FloatMatrix covariantMatrix;
  // double n = tangent_vectors.size();
  double n = surface_dimension;
  const FloatArray& rho_1 = tangent_vectors[0];
  if(n == 1) {
    covariantMatrix = {{rho_1.dotProduct(rho_1)}};
  } else if(n == 2) {
    const FloatArray& rho_2 = tangent_vectors[1];
    covariantMatrix = {
      {rho_1.dotProduct(rho_1), rho_1.dotProduct(rho_2)},
      {rho_2.dotProduct(rho_1), rho_2.dotProduct(rho_2)}
    };
  } else {
    OOFEM_WARNING("Incorect size of tangent_vectors");
  }
  return covariantMatrix;

}

FloatMatrix
StructuralPenaltyContactBoundaryCondition ::computeContravariantMetric(const std::vector<FloatArray> &tangent_vectors)
{
	FloatMatrix covariantMetric = this->computeCovariantMetric(tangent_vectors);
	FloatMatrix contravariantMetric;
	contravariantMetric.beInverseOf(covariantMetric);
	return contravariantMetric;
}

void
StructuralPenaltyContactBoundaryCondition :: setupContactSearchAlgorithm()
{
  if(this->surface_dimension == 2) {
    if (algo == 1) {
      using T = ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune;
      this->contactSearchAlgorithm = std::make_unique<T>(this->slaveContactSurface, this->masterContactSurface, domain);
    } else {
      using T = ContactSearchAlgorithm_Surface2FESurface_3d;
      this->contactSearchAlgorithm = std::make_unique<T>(this->slaveContactSurface, this->masterContactSurface, domain);
	}
  } else if(surface_dimension == 1) {
    this->contactSearchAlgorithm = std::make_unique<ContactSearchAlgorithm_Surface2FESurface_2d>(this->slaveContactSurface, this->masterContactSurface, domain);
  } else {
    OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition ::Unknown number of spatial dimensions");
  }
}


  
void
StructuralPenaltyContactBoundaryCondition::computeInternalForcesFromContact (FloatArray &answer, ContactPair *contactPair, TimeStep *tStep)
{
    answer.zero();
    auto normal_gap = contactPair->giveNormalGap();
    if ( normal_gap <= 1.e-15 ) {
      contactPair->initContactPoint();
      auto normal = contactPair->giveNormalVector();
      auto area = normal.computeNorm();
      normal /= area;
      FloatMatrix N;
      contactPair->computeNmatrix(N);
      // tractions in local coordinate system
      double normalTraction;
      FloatArray tangentialTraction, tangentialTractionTrial;
      ContactProcess mode = ContactProcess::None;
      this->computeTractions(normalTraction, tangentialTraction, tangentialTractionTrial, mode, contactPair, tStep);
      // traction in global coordinate system
      FloatArray tractionVector = normalTraction * normal;
      if (frictionShouldBeConsidered(friction, tStep)) {
		const auto tangentVectors = contactPair->giveTangentVectors();
		FloatMatrix contravariant_metric = this->computeContravariantMetric(tangentVectors);
		for (int i = 0; i < surface_dimension; i++){
			FloatArray rho_i = tangentVectors[i];
			for (int j = 0; j < surface_dimension; j++){
				double t_j = tangentialTraction(j);
				double a_ij = contravariant_metric(i, j);
				tractionVector += t_j * rho_i * a_ij;
			}
		}
      }
      answer.beTProductOf(N, tractionVector);
      contactPair->setTempTractionVector(tangentialTraction);
      //answer.times(area);
    } else {
      answer.clear();
    }
}


void
StructuralPenaltyContactBoundaryCondition :: computeTractions(double& normalTraction, FloatArray &tangentialTraction, FloatArray &tangentialTractionTrial, ContactProcess& mode, ContactPair* contactPair, TimeStep *tStep)
{
	// Algorithm from page 141
	// Computes normal and tangential tractions in local coordinate system
	tangentialTraction.zero();
	tangentialTraction.resize(surface_dimension);
	double normal_gap = contactPair->giveNormalGap();
	normalTraction = 0.0;
	if ( normal_gap > 0.0 ) {
		return;
	}
	// compute normal traction
	normalTraction = penalty_normal * normal_gap;
	if (frictionShouldBeConsidered(friction, tStep)) {
	  // auxiliary values
	  auto prevTangentialTractions = contactPair->giveTractionVector();
	  const auto prevTangentVectors = contactPair->givePreviousTangentVectors();
	  const auto tangentVectors = contactPair->giveTangentVectors();
	  FloatMatrix prev_contravariant_metric = this->computeContravariantMetric(prevTangentVectors);
	  FloatMatrix contravariant_metric = this->computeContravariantMetric(tangentVectors);
	  FloatMatrix covariant_metric = this->computeCovariantMetric(tangentVectors);
	  //
	  FloatArray delta_rho = contactPair->computeContactPointDisplacement();
	  //
	  for(int i = 0; i < surface_dimension; i++) {
		const auto& rho_i = tangentVectors[i];
	    for(int j = 0; j < surface_dimension; j++) {
	      const double a_ij = covariant_metric(i, j);
		  const auto& prev_rho_j = prevTangentVectors[j];
		  double prev_rho_j_dot_rho_i = prev_rho_j.dotProduct(rho_i);
	      for(int k = 0; k < surface_dimension; k++) {
		    const double prev_a_kj = prev_contravariant_metric(k, j);
		    const double prev_t_k = prevTangentialTractions(k);
		    tangentialTraction(i) += prev_t_k * prev_a_kj * prev_rho_j_dot_rho_i;
	      }
	      double delta_xi_j = 0.0;
	      for (int k = 0; k < surface_dimension; k++) {
		    const double a_jk = contravariant_metric(j, k);
		    const FloatArray& rho_k = tangentVectors[k];
		    delta_xi_j += delta_rho.dotProduct(rho_k) * a_jk;
	      }
	      tangentialTraction(i) -= penalty_tangential * delta_xi_j * a_ij;
	    }
	  }
	  // return mapping
	  double tNorm2 = 0.0;
	  for(int i = 0; i < surface_dimension; i++) {
	    double t_i = tangentialTraction(i);
	    for(int j = 0; j < surface_dimension; j++) {
	      double a_ij = contravariant_metric(i, j);
	      double t_j = tangentialTraction(j);
	      tNorm2 += t_i * t_j * a_ij;
	    }
	  }
	  double tNorm = std::sqrt(tNorm2);
	  double absNormalTraction = std::abs(normalTraction);
	  double f = tNorm - friction * absNormalTraction;
	  tangentialTractionTrial = tangentialTraction;
	  if (f > 0) {
		tangentialTractionTrial.times(-1);
		tangentialTraction.times(-1);
		double factor = friction * absNormalTraction / tNorm;
	    tangentialTraction *= factor;
	    mode = ContactProcess::Sliding;
	  } else {
	    mode = ContactProcess::Sticking;
	  }
	}
}



void
StructuralPenaltyContactBoundaryCondition::initializeFrom(InputRecord &ir)
{
    ContactBoundaryCondition::initializeFrom(ir);
    // contact surfaces
    IR_GIVE_FIELD(ir, this->slaveSurfaceNumber, _IFT_StructuralPenaltyContactBoundaryCondition_slaveSurfaceNum);
    IR_GIVE_FIELD(ir,this->masterSurfaceNumber, _IFT_StructuralPenaltyContactBoundaryCondition_masterSurfaceNum);
    //
    this->slaveContactSurface = dynamic_cast<StructuralFEContactSurface*>( this->giveDomain()->giveContactSurface(slaveSurfaceNumber));
    this->masterContactSurface = dynamic_cast<StructuralFEContactSurface*>( this->giveDomain()->giveContactSurface(masterSurfaceNumber));
    //
    IR_GIVE_FIELD(ir, this->penalty_normal, _IFT_StructuralPenaltyContactBoundaryCondition_penaltyNormal);
    IR_GIVE_FIELD(ir, this->penalty_tangential, _IFT_StructuralPenaltyContactBoundaryCondition_penaltyTangential);
    IR_GIVE_FIELD(ir, this->friction, _IFT_StructuralPenaltyContactBoundaryCondition_friction);
    int nsd;
    IR_GIVE_OPTIONAL_FIELD(ir, nsd, _IFT_StructuralPenaltyContactBoundaryCondition_nsd);
    surface_dimension = nsd - 1;
    //
    int algo = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, algo, _IFT_StructuralPenaltyContactBoundaryCondition_algo);
    this->algo = algo;
}


void
StructuralPenaltyContactBoundaryCondition :: postInitialize()
{
  // take nodes of node set after everything is initialized
  this->masterSurfaceElements = domain->giveSet(this->masterSurfaceNumber)->giveElementList();
  this->slaveSurfaceElements = domain->giveSet(this->slaveSurfaceNumber)->giveElementList();
  ContactBoundaryCondition :: postInitialize();
}



void
StructuralPenaltyContactBoundaryCondition :: giveLocationArrays(std::vector< IntArray > &rows, std::vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    //returns all possible combinations of dof that can theoretically be triggered by contact
    //of any segment with any node. There is a plenty room for optimization of this
    IntArray m_loc, s_loc;

    rows.resize(0);
    cols.resize(0);

    const auto& contactPairs = getContactPairs();
    for(auto const &cp : contactPairs) {
      if(cp->inContact()) {
	cp->giveSlaveContactPoint()->giveLocationArray( s_loc, dofs, c_s);
	cp->giveMasterContactPoint()->giveLocationArray( m_loc, dofs, r_s);
	// insert location arrays into the answer arrays
	rows.push_back(s_loc);
	cols.push_back(m_loc);
	rows.push_back(m_loc);
	cols.push_back(s_loc);
      }
    }
}


  

} // namespace oofem
