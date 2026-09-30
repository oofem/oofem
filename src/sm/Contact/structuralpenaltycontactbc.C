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
#include "floatarrayf.h"
#include "floatmatrix.h"
#include "classfactory.h"
#include "mathfem.h"
#include "sm/Elements/structuralelement.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>
#include "Contact/contactsearch.h"
#include "Contact/contactsearchsweepandprune.h"
#include "Contact/contactelement.h"
#include "sm/Materials/structuralmaterial.h"
#include "node.h"
#include "dof.h"
#include "engngm.h"
#include "timestep.h"
#include "datastream.h"
#include "contextioerr.h"
namespace oofem {
REGISTER_BoundaryCondition(StructuralPenaltyContactBoundaryCondition);

namespace {

void writeMatrixBlock(std::ostream &stream, const char *label, const FloatMatrix &matrix)
{
  stream << label << '\n';
  stream << matrix.giveNumberOfRows() << ' ' << matrix.giveNumberOfColumns() << '\n';
  stream << std::setprecision(16);
  for (int i = 1; i <= matrix.giveNumberOfRows(); ++i) {
    for (int j = 1; j <= matrix.giveNumberOfColumns(); ++j) {
      if (j > 1) {
        stream << ' ';
      }
      stream << matrix.at(i, j);
    }
    stream << '\n';
  }
}

void writeVectorBlock(std::ostream &stream, const char *label, const FloatArray &vector)
{
  stream << label << '\n';
  stream << vector.giveSize() << '\n';
  stream << std::setprecision(16);
  for (int i = 1; i <= vector.giveSize(); ++i) {
    if (i > 1) {
      stream << ' ';
    }
    stream << vector.at(i);
  }
  stream << '\n';
}

}


/** temporary auxiliary function */
bool frictionShouldBeConsidered(double friction, const ContactPair *contactPair) {
	// The discrete return map advances friction from the last converged active
	// state.  A projection retained while the pair was open is search history,
	// not a tangential-traction history: using it here would accumulate relative
	// motion before contact and create friction on the first closing increment.
	return friction > 0 && contactPair->hasContactHistory();
}



void StructuralPenaltyContactBoundaryCondition::computeTangentFromContact(FloatMatrix &answer, ContactPair *contactPair, TimeStep *tStep)
{
    answer.zero();
    if (contactPair->hasActiveContact()) {
      // The convected friction residual after a facet switch depends on the
      // current coordinates of both the new and committed master facets.  The
      // exact discrete tangent below supplies those rectangular history
      // columns for surface projections. Generalized edge/vertex and
      // directional-projection branches retain the numerical fallback.
      const bool generalizedFeatureTangent =
        generalizedFeatures
        && contactPair->giveCurrentMasterFeatureType() != ContactFeatureType::Surface;
      const bool frictionalFacetTransition =
        frictionShouldBeConsidered(friction, contactPair)
        && contactPair->hasMasterFeatureTransition();
      // The closed-form matrix below is the covariant rate-form tangent.  It
      // does not yet contain the complete algorithmic derivative of the
      // finite-step friction update from the committed material projection.
      // The automatic path uses the exact finite-step derivative below on a
      // smooth surface branch (including facet-history columns) and retains
      // the branch-frozen numerical derivative for unsupported projection
      // features. tangentmode 1 remains available for comparison with the
      // analytical rate form.
      const bool historyDependentFrictionTangent =
        frictionShouldBeConsidered(friction, contactPair);
      const bool automaticFiniteDifferenceTangent =
        directionalProjection || generalizedFeatureTangent
        || frictionalFacetTransition || historyDependentFrictionTangent;
      const bool smoothSurfaceDiscreteTangent =
        historyDependentFrictionTangent
        && !directionalProjection
        && !generalizedFeatureTangent
        && !contactPair->hasMasterFacetTransition()
        && contactPair->giveCurrentMasterFeatureType()
             == ContactFeatureType::Surface
        && contactPair->givePreviousMasterFeatureType()
             == ContactFeatureType::Surface;

      // tangentmode 3 explicitly requests the same exact derivative used by
      // the automatic path.  Its present scope is deliberately limited to a
      // smooth surface projection branch.  The discrete matrix includes
      // rectangular old-master columns after a facet transition; unsupported
      // generalized-feature and directional-projection cases retain the
      // branch-frozen numerical Jacobian.
      if ((tangentMode == 0 || tangentMode == 3)
          && smoothSurfaceDiscreteTangent
          && !finiteDifferenceCheckInProgress) {
        contactPair->initContactPoint();
        if (this->computeDiscreteTangent(answer, contactPair, tStep)) {
          if (finiteDifferenceCheck) {
            this->dumpFiniteDifferenceCheck(answer, contactPair, tStep);
          }
          return;
        }
      }
      const bool useFiniteDifferenceTangent =
        tangentMode == 2
        || (tangentMode == 3 && historyDependentFrictionTangent)
        || (tangentMode == 0 && automaticFiniteDifferenceTangent);
      FloatMatrix finiteDifferenceSolverTangent;
      bool returnFiniteDifferenceSolverTangent = false;
      if (useFiniteDifferenceTangent && !finiteDifferenceCheckInProgress) {
        finiteDifferenceSolverTangent = this->computeFiniteDifferenceTangent(
          this->findContactPairIndex(contactPair), tStep);
        const bool canAuditSmoothLocalRateForm =
          finiteDifferenceCheck
          && !directionalProjection
          && !generalizedFeatureTangent
          && !contactPair->hasMasterFacetTransition()
          && contactPair->giveResidualNodes()
               == contactPair->giveLinearizationNodes();
        if (!canAuditSmoothLocalRateForm) {
          answer = finiteDifferenceSolverTangent;
          return;
        }
        // Keep the nonlinear trajectory on the requested branch-frozen FD
        // Jacobian while still evaluating the closed-form matrix below for
        // the component audit.  Without this diagnostic path tangentmode 2
        // returned before fdcheck could report any quantities.
        returnFiniteDifferenceSolverTangent = true;
      }
      contactPair->initContactPoint();
      auto normal = contactPair->giveNormalVector();
      auto normalNorm = normal.computeNorm();
      if (normalNorm <= 0.0) {
        OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: zero contact normal");
      }
      normal /= normalNorm;
      // get the contact element shape function values at this node
      FloatMatrix N;
      contactPair->computeNmatrix(N);
      //
      ContactTractionState tractionState = this->computeTractionState(contactPair, tStep);
      double normal_traction = tractionState.normalTraction;
      FloatArray tangential_traction = tractionState.tangentialTraction;
      FloatArray tangential_traction_trial = tractionState.trialTangentialTraction;
      ContactProcess mode = tractionState.mode;
	  const double normalPenalty = this->giveNormalPenalty(contactPair);
	  const double tangentialPenalty = this->giveTangentialPenalty(contactPair);
	  double abs_normal_traction = std::abs(normal_traction);
      /****************** Normal part of the stiffness matrix ******************************/
      //@todo: update
      double dA = this->giveContactAreaForLinearization(contactPair, tStep);
      // Equation(7.17) page 189
      answer.plus_Nt_a_otimes_b_B(N, normal, normal, N, normalPenalty * dA);
      //
      auto tangent_vectors = contactPair->giveTangentVectors();
      auto contravariant_metric = this->computeContravariantMetric(tangent_vectors);
      // calculate curvature tensor kappa
      FloatMatrix curvatureCovariant;
      contactPair->computeCurvature(curvatureCovariant, tStep);
      FloatMatrix curvatureContravariant(surface_dimension, surface_dimension);
      curvatureContravariant.zero();
      for (int i = 0; i < surface_dimension; ++i) {
        for (int j = 0; j < surface_dimension; ++j) {
          for (int k = 0; k < surface_dimension; ++k) {
            for (int l = 0; l < surface_dimension; ++l) {
              curvatureContravariant(i, j) += contravariant_metric(i, k)
                                              * curvatureCovariant(k, l)
                                              * contravariant_metric(l, j);
            }
          }
        }
      }
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
		  double h_ij = curvatureContravariant(i, j);
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
      if (frictionShouldBeConsidered(friction, contactPair)) {
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
			  double h_ij = curvatureContravariant(i, j);
	          // Equation (7.26), page 192, with the sign matching OOFEM's assembled residual convention.
              factor = +tangentialPenalty * a_ij;
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
			  const double h_ij = curvatureContravariant(i,j);
			  // Equation (7.35a), page 194
			  factor = -normalPenalty * friction * t_i / t_norm * a_ij;
			  k_sl_constitutive_non_symmetric.plus_Nt_a_otimes_b_B(N, rho_j, normal, N, factor);
			  // Equation (7.35b), page 194
			  factor = -tangentialPenalty * friction * abs_normal_traction / t_norm * a_ij;
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
			      factor = +tangentialPenalty * friction * abs_normal_traction * t_i * t_j / t_norm3 * a_ik * a_jl;
			      k_sl_constitutive_symmetric_2.plus_Nt_a_otimes_b_B(N, rho_k, rho_l, N, factor);
			      // Equation (7.35d), page 194
			      const double factorCommon = -friction * abs_normal_traction * t_i / t_norm;
			      factor = factorCommon * a_il * a_jk;
			      k_sl_rotational_symmetric.plus_Nt_a_otimes_b_B(N, rho_k, rho_l, dNj, factor);
			      factor = factorCommon * a_ik * a_jl;
			      k_sl_rotational_symmetric.plus_Nt_a_otimes_b_B(dNj, rho_k, rho_l, N, factor);
			    }
	    	  }
	        }
	      }
              if (surface_dimension == 2) {
                std::vector<std::vector<FloatArray>> secondBaseVectors;
                contactPair->computeSecondBaseVectors(secondBaseVectors, tStep);

                std::vector<FloatArray> contravariantTangents(surface_dimension);
                FloatArray trialTractionContravariant(surface_dimension);
                trialTractionContravariant.zero();
                FloatArray globalTrialTraction(normal.giveSize());
                globalTrialTraction.zero();
                for (int i = 0; i < surface_dimension; ++i) {
                  contravariantTangents[i].resize(normal.giveSize());
                  contravariantTangents[i].zero();
                  for (int j = 0; j < surface_dimension; ++j) {
                    contravariantTangents[i].add(contravariant_metric(i, j), tangent_vectors[j]);
                    trialTractionContravariant(i) += contravariant_metric(i, j) * tangential_traction(j);
                  }
                  globalTrialTraction.add(tangential_traction(i), contravariantTangents[i]);
                }

                std::vector<double> connection(surface_dimension, 0.0);
                double curvatureContraction = 0.0;
                for (int j = 0; j < surface_dimension; ++j) {
                  for (int n = 0; n < surface_dimension; ++n) {
                    connection[n] += trialTractionContravariant(j)
                                     * globalTrialTraction.dotProduct(secondBaseVectors[j][n]);
                  }
                  for (int p = 0; p < surface_dimension; ++p) {
                    curvatureContraction += trialTractionContravariant(j)
                                            * trialTractionContravariant(p)
                                            * curvatureCovariant(p, j);
                  }
                }

                FloatArray rightVector = normal;
                rightVector.times(curvatureContraction);
                for (int n = 0; n < surface_dimension; ++n) {
                  rightVector.add(-connection[n], contravariantTangents[n]);
                }

                k_sl_curvature_non_symmetric.plus_Nt_a_otimes_b_B(
                  N, globalTrialTraction, rightVector, N,
                  friction * abs_normal_traction / t_norm3);
              } else {
                const double a11 = tangent_vectors[0].dotProduct(tangent_vectors[0]);
                const double oneDimensionalCurvatureFactor =
                  friction * abs_normal_traction
                  * std::copysign(1.0, tangential_traction(0))
                  * curvatureCovariant(0, 0) / std::pow(a11, 1.5);
                k_sl_curvature_non_symmetric.plus_Nt_a_otimes_b_B(
                  N, tangent_vectors[0], normal, N, oneDimensionalCurvatureFactor);
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
      if (finiteDifferenceCheck && !finiteDifferenceCheckInProgress) {
        this->dumpFiniteDifferenceCheck(answer, contactPair, tStep);
      }
      if (returnFiniteDifferenceSolverTangent) {
        answer = finiteDifferenceSolverTangent;
      }
    }
}

  
FloatMatrix
StructuralPenaltyContactBoundaryCondition ::computeCovariantMetric(const std::vector<FloatArray> &tangent_vectors)
{
  FloatMatrix covariantMatrix;
  if (static_cast<int>(tangent_vectors.size()) != surface_dimension) {
    OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: tangent vector count does not match contact surface dimension");
  }
  const int n = surface_dimension;
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

FloatArray
StructuralPenaltyContactBoundaryCondition::computeTransportedPreviousPhysicalTraction(
    ContactPair *contactPair)
{
  const int spatialDimension = surface_dimension + 1;
  const std::vector<FloatArray> previousTangents =
    contactPair->givePreviousTangentVectors();
  const FloatMatrix previousContravariant =
    this->computeContravariantMetric(previousTangents);
  const FloatArray previousTraction = contactPair->giveTractionVector();

  FloatArray transported(spatialDimension);
  transported.zero();
  for (int i = 0; i < surface_dimension; ++i) {
    for (int j = 0; j < surface_dimension; ++j) {
      transported.add(
        previousContravariant.at(i + 1, j + 1)
          * previousTraction.at(j + 1),
        previousTangents[i]);
    }
  }

  // Across a facet or feature switch there is no common smooth parametric
  // basis. Directly projecting the old vector onto the new tangent plane
  // loses magnitude and depends on the arbitrary edge/vertex frame. Transport
  // it with the unique minimal rotation that maps the committed normal to the
  // current normal. Transition branches use the numerical tangent, so the
  // derivative of this rotation is included by branch-frozen differences.
  if (!contactPair->hasMasterFeatureTransition()) {
    return transported;
  }

  FloatArray previousNormal = contactPair->givePreviousNormalVector();
  FloatArray currentNormal = contactPair->giveNormalVector();
  const double previousNorm = previousNormal.computeNorm();
  const double currentNorm = currentNormal.computeNorm();
  if (previousNormal.giveSize() != spatialDimension
      || currentNormal.giveSize() != spatialDimension
      || !(previousNorm > 0.0) || !(currentNorm > 0.0)) {
    return transported;
  }
  previousNormal /= previousNorm;
  currentNormal /= currentNorm;
  const double cosine = std::clamp(
    previousNormal.dotProduct(currentNormal), -1.0, 1.0);
  constexpr double rotationTolerance = 64.0
    * std::numeric_limits<double>::epsilon();

  FloatArray rotated(spatialDimension);
  rotated.zero();
  if (spatialDimension == 2) {
    const double sine = previousNormal.at(1) * currentNormal.at(2)
      - previousNormal.at(2) * currentNormal.at(1);
    rotated.at(1) = cosine * transported.at(1) - sine * transported.at(2);
    rotated.at(2) = sine * transported.at(1) + cosine * transported.at(2);
  } else if (spatialDimension == 3) {
    FloatArray axis(3);
    axis.at(1) = previousNormal.at(2) * currentNormal.at(3)
      - previousNormal.at(3) * currentNormal.at(2);
    axis.at(2) = previousNormal.at(3) * currentNormal.at(1)
      - previousNormal.at(1) * currentNormal.at(3);
    axis.at(3) = previousNormal.at(1) * currentNormal.at(2)
      - previousNormal.at(2) * currentNormal.at(1);
    double sine = axis.computeNorm();
    if (sine > rotationTolerance) {
      axis /= sine;
    } else if (cosine < 0.0) {
      // Antiparallel normals do not define a unique minimal axis. Use the
      // committed first tangent so the choice is material-history based and
      // independent of the arbitrary current vertex frame.
      axis = previousTangents.front();
      axis.add(-axis.dotProduct(previousNormal), previousNormal);
      double axisNorm = axis.computeNorm();
      if (!(axisNorm > rotationTolerance)) {
        axis = std::abs(previousNormal.at(1)) < 0.9
          ? Vec3(1.0, 0.0, 0.0)
          : Vec3(0.0, 1.0, 0.0);
        axis.add(-axis.dotProduct(previousNormal), previousNormal);
        axisNorm = axis.computeNorm();
      }
      axis /= axisNorm;
      sine = 0.0;
    } else {
      return transported;
    }

    FloatArray axisCrossTraction(3);
    axisCrossTraction.at(1) = axis.at(2) * transported.at(3)
      - axis.at(3) * transported.at(2);
    axisCrossTraction.at(2) = axis.at(3) * transported.at(1)
      - axis.at(1) * transported.at(3);
    axisCrossTraction.at(3) = axis.at(1) * transported.at(2)
      - axis.at(2) * transported.at(1);
    rotated = cosine * transported;
    rotated.add(sine, axisCrossTraction);
    rotated.add((1.0 - cosine) * axis.dotProduct(transported), axis);
  } else {
    return transported;
  }

  // Remove roundoff in the current normal direction. Since the rotation maps
  // tangent planes, this is normally at machine precision.
  rotated.add(-rotated.dotProduct(currentNormal), currentNormal);
  return rotated;
}

bool
StructuralPenaltyContactBoundaryCondition::computeDiscreteTangent(
    FloatMatrix &answer, ContactPair *contactPair, TimeStep *tStep)
{
  // Exact derivative of the finite-step residual actually evaluated by
  // computeTractionState.  The previous traction, basis and material
  // projection are committed history and therefore remain fixed during this
  // Newton linearization.
  auto *masterPoint =
    dynamic_cast<FEContactPoint *>(contactPair->giveMasterContactPoint());
  ContactElement *previousMaster = contactPair->givePreviousMasterContactElement();
  if (masterPoint == nullptr || previousMaster == nullptr) {
    return false;
  }

  const int spatialDimension = surface_dimension + 1;
  FloatArray normal = contactPair->giveNormalVector();
  const double normalNorm = normal.computeNorm();
  if (normalNorm <= 0.0) {
    OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: zero contact normal in discrete tangent");
  }
  normal /= normalNorm;
  const std::vector<FloatArray> tangents = contactPair->giveTangentVectors();
  const FloatMatrix covariantMetric = this->computeCovariantMetric(tangents);
  const FloatMatrix contravariantMetric =
    this->computeContravariantMetric(tangents);

  FloatMatrix N;
  contactPair->computeNmatrix(N);
  const std::vector<Node *> residualNodes = contactPair->giveResidualNodes();
  const std::vector<Node *> linearizationNodes =
    contactPair->giveLinearizationNodes();
  const int nResidualColumns = N.giveNumberOfColumns();
  const int nColumns =
    static_cast<int>(linearizationNodes.size()) * spatialDimension;
  if (N.giveNumberOfRows() != spatialDimension
      || dofs.giveSize() != spatialDimension
      || nResidualColumns
           != static_cast<int>(residualNodes.size()) * spatialDimension) {
    return false;
  }
  FloatMatrix geometryN(spatialDimension, nColumns);
  geometryN.zero();
  geometryN.setSubMatrix(N, 1, 1);

  std::vector<FloatMatrix> dNs;
  contactPair->compute_dNdxi_matrices(dNs);
  if (static_cast<int>(dNs.size()) != surface_dimension) {
    return false;
  }
  std::vector<FloatMatrix> geometryDerivatives;
  geometryDerivatives.reserve(surface_dimension);
  for (const FloatMatrix &dN : dNs) {
    if (dN.giveNumberOfRows() != spatialDimension
        || dN.giveNumberOfColumns() != nResidualColumns) {
      return false;
    }
    FloatMatrix expanded(spatialDimension, nColumns);
    expanded.zero();
    expanded.setSubMatrix(dN, 1, 1);
    geometryDerivatives.emplace_back(std::move(expanded));
  }

  FloatMatrix curvature;
  contactPair->computeCurvature(curvature, tStep);
  std::vector<std::vector<FloatArray>> secondBaseVectors;
  contactPair->computeSecondBaseVectors(secondBaseVectors, tStep);
  FloatMatrix projectionHessian(surface_dimension, surface_dimension);
  for (int i = 1; i <= surface_dimension; ++i) {
    for (int j = 1; j <= surface_dimension; ++j) {
      projectionHessian.at(i, j) = covariantMetric.at(i, j)
        - contactPair->giveNormalGap() * curvature.at(i, j);
    }
  }
  FloatMatrix inverseProjectionHessian;
  inverseProjectionHessian.beInverseOf(projectionHessian);

  // Differentiate the closest-point stationarity equations to obtain the
  // current master coordinates with respect to every pair displacement DOF.
  FloatMatrix localRightHandSide(surface_dimension, nColumns);
  for (int i = 0; i < surface_dimension; ++i) {
    for (int col = 1; col <= nColumns; ++col) {
      double tangentTimesN = 0.0;
      double normalTimesDerivativeN = 0.0;
      for (int r = 1; r <= spatialDimension; ++r) {
        tangentTimesN += tangents[i].at(r) * geometryN.at(r, col);
        normalTimesDerivativeN += normal.at(r)
          * geometryDerivatives[i].at(r, col);
      }
      localRightHandSide.at(i + 1, col) = tangentTimesN
        - contactPair->giveNormalGap() * normalTimesDerivativeN;
    }
  }
  FloatMatrix localDerivative;
  localDerivative.beProductOf(inverseProjectionHessian,
                              localRightHandSide);

  std::vector<FloatMatrix> tangentDerivatives(
    surface_dimension, FloatMatrix(spatialDimension, nColumns));
  for (int i = 0; i < surface_dimension; ++i) {
    for (int r = 1; r <= spatialDimension; ++r) {
      for (int col = 1; col <= nColumns; ++col) {
        double value = -geometryDerivatives[i].at(r, col);
        for (int j = 0; j < surface_dimension; ++j) {
          value += secondBaseVectors[i][j].at(r)
            * localDerivative.at(j + 1, col);
        }
        tangentDerivatives[i].at(r, col) = value;
      }
    }
  }

  std::vector<FloatArray> contravariantTangents(surface_dimension);
  for (int i = 0; i < surface_dimension; ++i) {
    contravariantTangents[i].resize(spatialDimension);
    contravariantTangents[i].zero();
    for (int j = 0; j < surface_dimension; ++j) {
      contravariantTangents[i].add(
        contravariantMetric.at(i + 1, j + 1), tangents[j]);
    }
  }
  FloatMatrix normalDerivative(spatialDimension, nColumns);
  normalDerivative.zero();
  for (int col = 1; col <= nColumns; ++col) {
    for (int i = 0; i < surface_dimension; ++i) {
      double normalTimesTangentDerivative = 0.0;
      for (int r = 1; r <= spatialDimension; ++r) {
        normalTimesTangentDerivative += normal.at(r)
          * tangentDerivatives[i].at(r, col);
      }
      for (int r = 1; r <= spatialDimension; ++r) {
        normalDerivative.at(r, col) -= contravariantTangents[i].at(r)
          * normalTimesTangentDerivative;
      }
    }
  }

  FloatMatrix covariantDerivative(
    surface_dimension * surface_dimension, nColumns);
  FloatMatrix contravariantDerivative(
    surface_dimension * surface_dimension, nColumns);
  for (int i = 0; i < surface_dimension; ++i) {
    for (int j = 0; j < surface_dimension; ++j) {
      const int row = i * surface_dimension + j + 1;
      for (int col = 1; col <= nColumns; ++col) {
        double value = 0.0;
        for (int r = 1; r <= spatialDimension; ++r) {
          value += tangentDerivatives[i].at(r, col) * tangents[j].at(r)
            + tangents[i].at(r) * tangentDerivatives[j].at(r, col);
        }
        covariantDerivative.at(row, col) = value;
      }
    }
  }
  for (int i = 0; i < surface_dimension; ++i) {
    for (int j = 0; j < surface_dimension; ++j) {
      const int row = i * surface_dimension + j + 1;
      for (int col = 1; col <= nColumns; ++col) {
        double value = 0.0;
        for (int k = 0; k < surface_dimension; ++k) {
          for (int l = 0; l < surface_dimension; ++l) {
            value -= contravariantMetric.at(i + 1, k + 1)
              * covariantDerivative.at(k * surface_dimension + l + 1, col)
              * contravariantMetric.at(l + 1, j + 1);
          }
        }
        contravariantDerivative.at(row, col) = value;
      }
    }
  }

  // The slip chord joins the current material projection to its committed
  // material position evaluated in the current configuration.  Current and
  // committed master facets can differ, so scatter both interpolation
  // matrices by node identity into the full (possibly rectangular) column
  // layout. Shared-node contributions then add in the ordinary way.
  FloatMatrix currentMasterN;
  masterPoint->computeNmatrix(currentMasterN);
  FloatMatrix previousMasterN;
  previousMaster->computeNmatrixAt(
    contactPair->givePreviousMasterLocalCoordinates(), previousMasterN);
  ContactElement *currentMaster = masterPoint->giveContactElement();
  if (currentMaster == nullptr
      || currentMasterN.giveNumberOfRows() != spatialDimension
      || previousMasterN.giveNumberOfRows() != spatialDimension
      || currentMasterN.giveNumberOfColumns()
           != currentMaster->giveNumberOfNodes() * spatialDimension
      || previousMasterN.giveNumberOfColumns()
           != previousMaster->giveNumberOfNodes() * spatialDimension) {
    return false;
  }
  FloatMatrix chordDerivative(spatialDimension, nColumns);
  chordDerivative.zero();

  const auto scatterMasterInterpolation =
    [&](ContactElement *element, const FloatMatrix &interpolation,
        double scale) {
      for (int nodeIndex = 1; nodeIndex <= element->giveNumberOfNodes(); ++nodeIndex) {
        Node *node = element->giveNode(nodeIndex);
        const auto position =
          std::find(linearizationNodes.begin(), linearizationNodes.end(), node);
        if (position == linearizationNodes.end()) {
          return false;
        }
        const int destinationOffset =
          static_cast<int>(position - linearizationNodes.begin())
          * spatialDimension;
        const int sourceOffset = (nodeIndex - 1) * spatialDimension;
        for (int r = 1; r <= spatialDimension; ++r) {
          for (int component = 1; component <= spatialDimension; ++component) {
            chordDerivative.at(r, destinationOffset + component) += scale
              * interpolation.at(r, sourceOffset + component);
          }
        }
      }
      return true;
    };
  if (!scatterMasterInterpolation(currentMaster, currentMasterN, 1.0)
      || !scatterMasterInterpolation(previousMaster, previousMasterN, -1.0)) {
    return false;
  }

  for (int r = 1; r <= spatialDimension; ++r) {
    for (int col = 1; col <= nColumns; ++col) {
      for (int i = 0; i < surface_dimension; ++i) {
        chordDerivative.at(r, col) += tangents[i].at(r)
          * localDerivative.at(i + 1, col);
      }
    }
  }

  const ContactTractionState tractionState =
    this->computeTractionState(contactPair, tStep);
  const FloatArray transportedPreviousPhysical =
    this->computeTransportedPreviousPhysicalTraction(contactPair);
  const FloatArray projectionChord =
    contactPair->computeContactPointDisplacement(tStep);

  FloatMatrix trialDerivative(surface_dimension, nColumns);
  trialDerivative.zero();
  const double tangentialPenalty = this->giveTangentialPenalty(contactPair);
  for (int i = 0; i < surface_dimension; ++i) {
    for (int col = 1; col <= nColumns; ++col) {
      double transportedDerivative = 0.0;
      double projectedChordDerivative = 0.0;
      for (int r = 1; r <= spatialDimension; ++r) {
        transportedDerivative += transportedPreviousPhysical.at(r)
          * tangentDerivatives[i].at(r, col);
        projectedChordDerivative += chordDerivative.at(r, col)
            * tangents[i].at(r)
          + projectionChord.at(r) * tangentDerivatives[i].at(r, col);
      }
      trialDerivative.at(i + 1, col) = transportedDerivative
        + tangentialPenalty * projectedChordDerivative;
    }
  }

  FloatMatrix gapDerivative(1, nColumns);
  FloatMatrix normalTractionDerivative(1, nColumns);
  const double normalPenalty = this->giveNormalPenalty(contactPair);
  for (int col = 1; col <= nColumns; ++col) {
    double value = 0.0;
    for (int r = 1; r <= spatialDimension; ++r) {
      value += normal.at(r) * geometryN.at(r, col);
    }
    gapDerivative.at(1, col) = value;
    normalTractionDerivative.at(1, col) = normalPenalty * value;
  }

  FloatMatrix tangentialNormDerivative(1, nColumns);
  tangentialNormDerivative.zero();
  if (tractionState.tangentialNorm > 0.0) {
    for (int col = 1; col <= nColumns; ++col) {
      double squaredNormDerivative = 0.0;
      for (int i = 0; i < surface_dimension; ++i) {
        for (int j = 0; j < surface_dimension; ++j) {
          const double aij = contravariantMetric.at(i + 1, j + 1);
          const double trialI = tractionState.trialTangentialTraction.at(i + 1);
          const double trialJ = tractionState.trialTangentialTraction.at(j + 1);
          squaredNormDerivative += trialDerivative.at(i + 1, col) * aij * trialJ;
          squaredNormDerivative += trialI
            * contravariantDerivative.at(i * surface_dimension + j + 1, col)
            * trialJ;
          squaredNormDerivative += trialI * aij
            * trialDerivative.at(j + 1, col);
        }
      }
      tangentialNormDerivative.at(1, col) =
        0.5 * squaredNormDerivative / tractionState.tangentialNorm;
    }
  }

  FloatMatrix returnedDerivative(surface_dimension, nColumns);
  returnedDerivative.zero();
  FloatMatrix slidingDerivative(surface_dimension, nColumns);
  slidingDerivative.zero();
  if (tractionState.tangentialNorm > 0.0) {
    const double normalTractionSign =
      tractionState.normalTraction < 0.0 ? -1.0
      : (tractionState.normalTraction > 0.0 ? 1.0 : 0.0);
    const double returnScale =
      tractionState.frictionBound / tractionState.tangentialNorm;
    for (int i = 1; i <= surface_dimension; ++i) {
      for (int col = 1; col <= nColumns; ++col) {
        const double frictionBoundDerivative =
          friction * normalTractionSign
            * normalTractionDerivative.at(1, col);
        const double scaleDerivative =
          frictionBoundDerivative / tractionState.tangentialNorm
          - tractionState.frictionBound
            * tangentialNormDerivative.at(1, col)
            / (tractionState.tangentialNorm * tractionState.tangentialNorm);
        slidingDerivative.at(i, col) = returnScale
            * trialDerivative.at(i, col)
          + tractionState.trialTangentialTraction.at(i) * scaleDerivative;
      }
    }
  }

  // The sliding limit for isotropic hardening is the weighted sum of the
  // elastic trial derivative and the perfect-Coulomb return derivative.
  FloatMatrix limitingSlidingDerivative = slidingDerivative;
  if (frictionHardening > 0.0) {
    limitingSlidingDerivative.times(1.0 - frictionHardening);
    limitingSlidingDerivative.add(frictionHardening, trialDerivative);
  }

  if (tractionState.mode == ContactProcess::Sticking) {
    returnedDerivative = trialDerivative;
  } else if (tractionState.mode == ContactProcess::Sliding) {
    returnedDerivative = slidingDerivative;
  }

  // A differentiable return projection prevents groups of integration points
  // from alternating indefinitely between the two sharp Coulomb branches.
  // The regularized residual is radial in the trial traction, so its exact
  // algorithmic derivative follows from the returned traction magnitude.
  if ((frictionTransition > 0.0 || frictionHardening > 0.0)
      && tractionState.tangentialNorm > 0.0) {
    double returnedMagnitude = 0.0;
    double magnitudeDerivativeWrtTrial = 0.0;
    double magnitudeDerivativeWrtBound = 0.0;
    this->computeFrictionReturnMagnitude(
      returnedMagnitude, magnitudeDerivativeWrtTrial,
      magnitudeDerivativeWrtBound, tractionState.tangentialNorm,
      tractionState.frictionBound);
    const double returnedScale =
      returnedMagnitude / tractionState.tangentialNorm;
    for (int i = 1; i <= surface_dimension; ++i) {
      for (int col = 1; col <= nColumns; ++col) {
        const double frictionBoundDerivative =
          friction
          * (tractionState.normalTraction < 0.0 ? -1.0 : 1.0)
          * normalTractionDerivative.at(1, col);
        const double scaleDerivative =
          (magnitudeDerivativeWrtTrial * tractionState.tangentialNorm
             - returnedMagnitude)
            / (tractionState.tangentialNorm * tractionState.tangentialNorm)
            * tangentialNormDerivative.at(1, col)
          + magnitudeDerivativeWrtBound / tractionState.tangentialNorm
            * frictionBoundDerivative;
        returnedDerivative.at(i, col) =
          returnedScale * trialDerivative.at(i, col)
          + tractionState.trialTangentialTraction.at(i) * scaleDerivative;
      }
    }
  }

  // At phi = ||t_trial|| - mu |t_n| = 0, the Coulomb return map is
  // semismooth rather than classically differentiable.  Choosing either the
  // sticking or sliding limit for every point can make a whole contact patch
  // change branch in lockstep at the first Newton iteration.  The midpoint of
  // the two limiting Jacobians is a member of the Clarke generalized
  // Jacobian and is also the derivative selected by a centered difference at
  // this kink.
  const double yieldValue =
    tractionState.tangentialNorm - tractionState.frictionBound;
  const double yieldTolerance =
    std::sqrt(std::numeric_limits<double>::epsilon())
    * std::max({tractionState.tangentialNorm,
                tractionState.frictionBound, 1.e-12});
  if (frictionTransition == 0.0
      && tractionState.tangentialNorm > yieldTolerance
      && tractionState.frictionBound > yieldTolerance
      && std::abs(yieldValue) <= yieldTolerance) {
    if (tractionState.mode == ContactProcess::Sticking) {
      returnedDerivative.add(limitingSlidingDerivative);
    } else {
      returnedDerivative.add(trialDerivative);
    }
    returnedDerivative.times(0.5);
  }

  FloatArray physicalReturned(spatialDimension);
  physicalReturned.zero();
  FloatMatrix physicalReturnedDerivative(spatialDimension, nColumns);
  physicalReturnedDerivative.zero();
  for (int i = 0; i < surface_dimension; ++i) {
    for (int j = 0; j < surface_dimension; ++j) {
      const double aij = contravariantMetric.at(i + 1, j + 1);
      const double returnedJ = tractionState.tangentialTraction.at(j + 1);
      physicalReturned.add(aij * returnedJ, tangents[i]);
      for (int r = 1; r <= spatialDimension; ++r) {
        for (int col = 1; col <= nColumns; ++col) {
          physicalReturnedDerivative.at(r, col) +=
            tangentDerivatives[i].at(r, col) * aij * returnedJ;
          physicalReturnedDerivative.at(r, col) += tangents[i].at(r)
            * contravariantDerivative.at(i * surface_dimension + j + 1, col)
            * returnedJ;
          physicalReturnedDerivative.at(r, col) += tangents[i].at(r) * aij
            * returnedDerivative.at(j + 1, col);
        }
      }
    }
  }

  const FloatArray totalTraction =
    tractionState.normalTraction * normal + physicalReturned;
  FloatMatrix totalTractionDerivative(spatialDimension, nColumns);
  totalTractionDerivative.zero();
  for (int r = 1; r <= spatialDimension; ++r) {
    for (int col = 1; col <= nColumns; ++col) {
      totalTractionDerivative.at(r, col) =
        normalTractionDerivative.at(1, col) * normal.at(r)
        + tractionState.normalTraction * normalDerivative.at(r, col)
        + physicalReturnedDerivative.at(r, col);
    }
  }

  const double contactArea =
    this->giveContactAreaForLinearization(contactPair, tStep);
  answer.resize(nResidualColumns, nColumns);
  answer.zero();
  for (int row = 1; row <= nResidualColumns; ++row) {
    for (int col = 1; col <= nColumns; ++col) {
      for (int r = 1; r <= spatialDimension; ++r) {
        double interpolationDerivative = 0.0;
        for (int i = 0; i < surface_dimension; ++i) {
          interpolationDerivative += dNs[i].at(r, row)
            * localDerivative.at(i + 1, col);
        }
        answer.at(row, col) += contactArea
          * (interpolationDerivative * totalTraction.at(r)
             + N.at(r, row) * totalTractionDerivative.at(r, col));
      }
    }
  }
  return true;
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

  auto surfaceSearch =
    dynamic_cast<ContactSearchAlgorithm_Surface2FESurface *>(contactSearchAlgorithm.get());
  if (!surfaceSearch) {
    OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: incompatible contact search algorithm");
  }
  surfaceSearch->setSearchPadding(searchPadding);
  surfaceSearch->setGeneralizedFeatures(generalizedFeatures);
  surfaceSearch->setDirectionalProjection(directionalProjection);
  surfaceSearch->setFacetOwnershipHysteresis(facetOwnershipHysteresis);
  masterContactSurface->setDomainTolerance(searchTol);
  slaveContactSurface->setDomainTolerance(searchTol);
}


  
void
StructuralPenaltyContactBoundaryCondition::computeInternalForcesFromContact (FloatArray &answer, ContactPair *contactPair, TimeStep *tStep)
{
    answer.zero();
    if (contactPair->hasActiveContact()) {
      contactPair->initContactPoint();
      auto normal = contactPair->giveNormalVector();
      auto normalNorm = normal.computeNorm();
      if (normalNorm <= 0.0) {
        OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: zero contact normal");
      }
      normal /= normalNorm;
      FloatMatrix N;
      contactPair->computeNmatrix(N);
      // tractions in local coordinate system
      ContactTractionState tractionState = this->computeTractionState(contactPair, tStep);
      double normalTraction = tractionState.normalTraction;
      FloatArray tangentialTraction = tractionState.tangentialTraction;
      if (!finiteDifferenceCheckInProgress) {
        int outputStatus = 1;
        if (tractionState.mode == ContactProcess::Sticking) {
          outputStatus = 2;
        } else if (tractionState.mode == ContactProcess::Sliding) {
          outputStatus = 3;
        }
        contactPair->setOutputContactState(-normalTraction, outputStatus);
      }
      // traction in global coordinate system
      FloatArray tractionVector = normalTraction * normal;
      if (frictionShouldBeConsidered(friction, contactPair)) {
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
      if (!finiteDifferenceCheckInProgress) {
        contactPair->setTempTractionVector(tangentialTraction);
        contactPair->setTempAccumulatedPlasticSlip(
          tractionState.accumulatedPlasticSlip);
      }
      answer.times(this->giveContactAreaForLinearization(contactPair, tStep));
    } else {
      answer.clear();
    }
}

double
StructuralPenaltyContactBoundaryCondition :: giveContactAreaForLinearization(ContactPair *contactPair, TimeStep *tStep) const
{
  if (finiteDifferenceCheckInProgress && contactPair == finiteDifferenceContactPair) {
    return finiteDifferenceContactArea;
  }
  return contactPair->giveContactArea(tStep);
}

double
StructuralPenaltyContactBoundaryCondition :: giveAutomaticPenaltyFactor(ContactPair *contactPair) const
{
  if (!autoPenalty) {
    return 1.0;
  }

  auto *slavePoint = dynamic_cast<FEContactPoint *>(contactPair->giveSlaveContactPoint());
  if (!slavePoint) {
    OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: automatic penalty requires an FE slave point");
  }
  const auto factor = automaticPenaltyFactors.find(slavePoint->giveContactElementId());
  if (factor == automaticPenaltyFactors.end()) {
    OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: automatic penalty is missing for slave contact element %d",
                slavePoint->giveContactElementId());
  }
  return factor->second;
}

double
StructuralPenaltyContactBoundaryCondition :: giveNormalPenalty(ContactPair *contactPair) const
{
  return penalty_normal * this->giveAutomaticPenaltyFactor(contactPair);
}

double
StructuralPenaltyContactBoundaryCondition :: giveTangentialPenalty(ContactPair *contactPair) const
{
  return penalty_tangential * this->giveAutomaticPenaltyFactor(contactPair);
}

void
StructuralPenaltyContactBoundaryCondition :: initializeAutomaticPenaltyFactors()
{
  automaticPenaltyFactors.clear();
  if (!autoPenalty) {
    return;
  }

  // The contact integration rule provides a reliable reference face area for
  // both triangles and quadrilaterals, including distorted faces.
  std::map<int, double> referenceFaceAreas;
  for (const auto &pair : getContactPairs()) {
    auto *slavePoint = dynamic_cast<FEContactPoint *>(pair->giveSlaveContactPoint());
    if (!slavePoint) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: automatic penalty requires FE contact pairs");
    }
    referenceFaceAreas[slavePoint->giveContactElementId()] += pair->giveContactArea(nullptr);
  }

  double minimumFactor = std::numeric_limits<double>::infinity();
  double maximumFactor = 0.0;
  for (int surfaceIndex = 1;
       surfaceIndex <= slaveContactSurface->giveNumberOfContactElements();
       ++surfaceIndex) {
    ContactElement *face = slaveContactSurface->giveContactElement_InSet(surfaceIndex);
    Element *parent = nullptr;

    // OOFEM contact faces are independent elements, so recover the adjacent
    // continuum element by the conforming face-node subset.  This is done once
    // during initialization, not during assembly.
    for (int elementIndex = 1; elementIndex <= domain->giveNumberOfElements(); ++elementIndex) {
      Element *candidate = domain->giveElement(elementIndex);
      if (dynamic_cast<ContactElement *>(candidate) != nullptr
          || candidate->giveNumberOfNodes() < face->giveNumberOfNodes()) {
        continue;
      }

      bool containsFace = true;
      for (int faceNode = 1; faceNode <= face->giveNumberOfNodes() && containsFace; ++faceNode) {
        bool foundNode = false;
        for (int candidateNode = 1; candidateNode <= candidate->giveNumberOfNodes(); ++candidateNode) {
          if (face->giveDofManager(faceNode) == candidate->giveDofManager(candidateNode)) {
            foundNode = true;
            break;
          }
        }
        containsFace = foundNode;
      }
      if (containsFace) {
        parent = candidate;
        break;
      }
    }

    if (!parent) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: no parent volume element found for contact element %d",
                  face->giveNumber());
    }

    // Elements whose material is assigned only through a Set-based cross
    // section (e.g. `SimpleCS ... material M set S`, common in generated
    // inputs) never populate Element::material, so giveMaterial() has
    // nothing to resolve and asserts/throws if called on them. Check the
    // raw material number first (never throws) and fall back to the
    // element's cross section, the same path StructuralElement material
    // access normally goes through.
    StructuralMaterial *material = parent->giveMaterialNumber() > 0
        ? dynamic_cast<StructuralMaterial *>(parent->giveMaterial()) : nullptr;
    if (!material) {
      if (auto *structuralParent = dynamic_cast<StructuralElement *>(parent)) {
        if (IntegrationRule *iRule = structuralParent->giveDefaultIntegrationRulePtr()) {
          if (GaussPoint *gp = iRule->getIntegrationPoint(0)) {
            if (StructuralCrossSection *cs = structuralParent->giveStructuralCrossSection()) {
              material = dynamic_cast<StructuralMaterial *>(cs->giveMaterial(gp));
            }
          }
        }
      }
    }
    if (!material) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: parent element %d has no resolvable structural material",
                  parent->giveNumber());
    }
    const double modulus = material->giveInitialYoungsModulus();
    // Use the element integration rule instead of Element::computeVolume();
    // several 3-D structural interpolations (including LSpace) deliberately
    // leave FEInterpolation3d::giveVolume unimplemented.
    const double volume = std::abs(parent->computeVolumeAreaOrLength());
    const auto faceArea = referenceFaceAreas.find(face->giveNumber());
    if (!(modulus > 0.0) || !(volume > 0.0)
        || faceArea == referenceFaceAreas.end() || !(faceArea->second > 0.0)) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: invalid automatic penalty data for contact element %d",
                  face->giveNumber());
    }

    const double factor = modulus * faceArea->second / volume;
    automaticPenaltyFactors[face->giveNumber()] = factor;
    minimumFactor = std::min(minimumFactor, factor);
    maximumFactor = std::max(maximumFactor, factor);
  }

  OOFEM_LOG_INFO("Contact BC %d: automatic penalty E_n*A/V range [%e, %e], pn/pt scales [%e, %e]\n",
                 this->giveNumber(), minimumFactor, maximumFactor,
                 penalty_normal, penalty_tangential);
}

double
StructuralPenaltyContactBoundaryCondition :: computeTangentialNorm(const FloatArray &tangentialTraction, const FloatMatrix &contravariantMetric) const
{
  double tNorm2 = 0.0;
  for (int i = 0; i < surface_dimension; i++) {
    const double t_i = tangentialTraction(i);
    for (int j = 0; j < surface_dimension; j++) {
      tNorm2 += t_i * tangentialTraction(j) * contravariantMetric(i, j);
    }
  }
  return std::sqrt(std::max(0.0, tNorm2));
}

void
StructuralPenaltyContactBoundaryCondition :: computeFrictionReturnMagnitude(
    double &returnedMagnitude,
    double &derivativeWrtTrialMagnitude,
    double &derivativeWrtFrictionBound,
    double trialMagnitude,
    double frictionBound) const
{
  returnedMagnitude = 0.0;
  derivativeWrtTrialMagnitude = 0.0;
  derivativeWrtFrictionBound = 0.0;
  if (trialMagnitude <= 0.0) {
    return;
  }

  double perfectlyPlasticMagnitude = 0.0;
  double perfectDerivativeWrtTrial = 0.0;
  double perfectDerivativeWrtBound = 1.0;
  if (frictionBound > 0.0) {
    const double ratio = trialMagnitude / frictionBound;
    if (frictionTransition == 0.0) {
      if (ratio <= 1.0) {
        perfectlyPlasticMagnitude = trialMagnitude;
        perfectDerivativeWrtTrial = 1.0;
        perfectDerivativeWrtBound = 0.0;
      } else {
        perfectlyPlasticMagnitude = frictionBound;
      }
    } else {
      // Smooth p-norm approximation of min(r,b):
      //   m = r / (1 + (r/b)^p)^(1/p),  p = 2/eta.
      // It approaches r in stick and b in slip without an active-set kink.
      const double exponent = 2.0 / frictionTransition;
      const double logarithmicRatio = exponent * std::log(ratio);
      double softplus = 0.0;
      double slidingWeight = 0.0;
      if (logarithmicRatio > 50.0) {
        softplus = logarithmicRatio + std::log1p(std::exp(-logarithmicRatio));
        slidingWeight = 1.0 / (1.0 + std::exp(-logarithmicRatio));
      } else if (logarithmicRatio < -50.0) {
        const double exponential = std::exp(logarithmicRatio);
        softplus = exponential;
        slidingWeight = exponential / (1.0 + exponential);
      } else {
        const double exponential = std::exp(logarithmicRatio);
        softplus = std::log1p(exponential);
        slidingWeight = exponential / (1.0 + exponential);
      }
      const double returnedScale = std::exp(-softplus / exponent);
      perfectlyPlasticMagnitude = trialMagnitude * returnedScale;
      perfectDerivativeWrtTrial = returnedScale * (1.0 - slidingWeight);
      perfectDerivativeWrtBound = ratio * returnedScale * slidingWeight;
    }
  }

  // Isotropic hardening with H = h*k_t/(1-h) gives the returned magnitude
  // h*r + (1-h)*min(r,b+H*alpha_n). This exposes h directly as the
  // post-yield tangent ratio and extends the smooth cap consistently.
  returnedMagnitude = frictionHardening * trialMagnitude
    + (1.0 - frictionHardening) * perfectlyPlasticMagnitude;
  derivativeWrtTrialMagnitude = frictionHardening
    + (1.0 - frictionHardening) * perfectDerivativeWrtTrial;
  derivativeWrtFrictionBound =
    (1.0 - frictionHardening) * perfectDerivativeWrtBound;
}

std::size_t
StructuralPenaltyContactBoundaryCondition::findContactPairIndex(const ContactPair *contactPair) const
{
  const auto &contactPairs = getContactPairs();
  for (std::size_t i = 0; i < contactPairs.size(); ++i) {
    if (contactPairs[i].get() == contactPair) {
      return i;
    }
  }

  OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: contact pair not found in active search list");
}

bool
StructuralPenaltyContactBoundaryCondition::updateFiniteDifferencePairGeometry(
    ContactPair *contactPair, int masterElementId,
    ContactFeatureType featureType, int featureIndex, TimeStep *tStep)
{
  ContactElement *masterElement = masterContactSurface->giveContactElement(masterElementId);
  ContactPoint *slavePoint = contactPair->giveSlaveContactPoint();

  if (surface_dimension == 2) {
    if (featureType != ContactFeatureType::Surface) {
      ContactProjection projection;
      if (featureType == ContactFeatureType::Edge) {
        projection = masterContactSurface->findContactPointOnEdge_3d(
          slavePoint, masterElement, featureIndex, tStep);
      } else {
        projection = masterContactSurface->findContactPointOnVertex_3d(
          slavePoint, masterElement, featureIndex, tStep);
      }
      if (!projection.valid) {
        return false;
      }
      auto masterPoint = std::make_unique<FEContactPoint_Master>(
        masterContactSurface, masterElementId, surface_dimension,
        projection.localCoordinates, featureType, featureIndex);
      contactPair->setMasterContactPoint(std::move(masterPoint));
      contactPair->setNormalGap(projection.gap);
      contactPair->setNormalVector(projection.normal);
      contactPair->setTangentVector1(projection.tangent1);
      contactPair->setTangentVector2(projection.tangent2);
      return true;
    }

    const FloatArrayF<3> direction = directionalProjection
      ? slaveContactSurface->computeContactPointNormal_3d(slavePoint, tStep)
      : FloatArrayF<3> {};
    if (directionalProjection && !(norm(direction) > 0.0)) {
      return false;
    }
    auto [inside, localCoordinates, gap, distanceSquared, normal, tangent1, tangent2] =
      directionalProjection
      ? masterContactSurface->findContactPointAlongDirectionInElement_3d(
          slavePoint, masterElement, direction, tStep)
      : masterContactSurface->findContactPointInElement_3d(
          slavePoint, masterElement, tStep);
    (void) distanceSquared;
    if (!inside) {
      return false;
    }
    auto masterPoint = std::make_unique<FEContactPoint_Master>(
      masterContactSurface, masterElementId, surface_dimension, localCoordinates);
    contactPair->setMasterContactPoint(std::move(masterPoint));
    contactPair->setNormalGap(gap);
    contactPair->setNormalVector(normal);
    contactPair->setTangentVector1(tangent1);
    contactPair->setTangentVector2(tangent2);
    return true;
  }

  if (featureType == ContactFeatureType::Vertex) {
    ContactProjection projection =
      masterContactSurface->findContactPointOnVertex_2d(
        slavePoint, masterElement, featureIndex, tStep);
    if (!projection.valid) {
      return false;
    }
    auto masterPoint = std::make_unique<FEContactPoint_Master>(
      masterContactSurface, masterElementId, surface_dimension,
      projection.localCoordinates, featureType, featureIndex);
    contactPair->setMasterContactPoint(std::move(masterPoint));
    contactPair->setNormalGap(projection.gap);
    contactPair->setNormalVector(projection.normal);
    contactPair->setTangentVector1(projection.tangent1);
    return true;
  }

  auto [inside, localCoordinates, gap, distanceSquared, normal, tangent1] =
    masterContactSurface->findContactPointInElement_2d(slavePoint, masterElement, tStep);
    (void) distanceSquared;
  if (!inside) {
    return false;
  }
  auto masterPoint = std::make_unique<FEContactPoint_Master>(
    masterContactSurface, masterElementId, surface_dimension, localCoordinates);
  contactPair->setMasterContactPoint(std::move(masterPoint));
  contactPair->setNormalGap(gap);
  contactPair->setNormalVector(normal);
  contactPair->setTangentVector1(tangent1);
  return true;
}

FloatMatrix
StructuralPenaltyContactBoundaryCondition::computeFiniteDifferenceTangent(
    std::size_t pairIndex, TimeStep *tStep)
{
  auto &contactPairs = getContactPairs();
  if (pairIndex >= contactPairs.size()) {
    OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: finite-difference pair index out of range");
  }

  ContactPair *contactPair = contactPairs[pairIndex].get();
  auto *masterPoint = dynamic_cast<FEContactPoint *>(contactPair->giveMasterContactPoint());
  if (masterPoint == nullptr) {
    OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: finite-difference check requires an FE master point");
  }
  const int masterElementId = masterPoint->giveContactElementId();
  const ContactFeatureType featureType = masterPoint->giveContactFeatureType();
  const int featureIndex = masterPoint->giveContactFeatureIndex();
  const std::vector<Node *> residualNodes = contactPair->giveResidualNodes();
  const std::vector<Node *> linearizationNodes = contactPair->giveLinearizationNodes();
  const int nRows = static_cast<int>(residualNodes.size()) * dofs.giveSize();
  const int nColumns = static_cast<int>(linearizationNodes.size()) * dofs.giveSize();
  FloatMatrix numericalTangent(nRows, nColumns);
  numericalTangent.zero();

  if (dofs.giveSize() != surface_dimension + 1) {
    OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: inconsistent finite-difference DOF layout");
  }

  std::vector<Dof *> perturbationDofs;
  perturbationDofs.reserve(nColumns);
  for (Node *node : linearizationNodes) {
    for (int component = 1; component <= dofs.giveSize(); ++component) {
      Dof *dof = node->giveDofWithID(dofs.at(component));
      if (dof == nullptr || !dof->isPrimaryDof()) {
        OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: finite-difference check requires primary displacement DOFs");
      }
      perturbationDofs.push_back(dof);
    }
  }

  std::vector<double> originalUnknowns;
  originalUnknowns.reserve(perturbationDofs.size());
  for (Dof *dof : perturbationDofs) {
    originalUnknowns.push_back(dof->giveUnknown(VM_Total, tStep));
  }

  finiteDifferenceContactPair = contactPair;
  finiteDifferenceContactArea = contactPair->giveContactArea(tStep);
  finiteDifferenceCheckInProgress = true;

  FloatArray baselineForce;
  this->computeInternalForcesFromContact(baselineForce, contactPair, tStep);
  baselineForce.resizeWithValues(nRows);

  for (int col = 1; col <= nColumns; ++col) {
    Dof *dof = perturbationDofs[col - 1];
    const double originalUnknown = originalUnknowns[col - 1];

    dof->updateUnknownsDictionary(tStep, VM_Total, originalUnknown + finiteDifferencePerturbation);
    const bool plusValid = this->updateFiniteDifferencePairGeometry(
      contactPair, masterElementId, featureType, featureIndex, tStep);
    FloatArray plusForce;
    if (plusValid) {
      this->computeInternalForcesFromContact(plusForce, contactPair, tStep);
      plusForce.resizeWithValues(nRows);
    }

    dof->updateUnknownsDictionary(tStep, VM_Total, originalUnknown - finiteDifferencePerturbation);
    const bool minusValid = this->updateFiniteDifferencePairGeometry(
      contactPair, masterElementId, featureType, featureIndex, tStep);
    FloatArray minusForce;
    if (minusValid) {
      this->computeInternalForcesFromContact(minusForce, contactPair, tStep);
      minusForce.resizeWithValues(nRows);
    }

    FloatArray column;
    if (plusValid && minusValid) {
      column.beDifferenceOf(plusForce, minusForce);
      column.times(0.5 / finiteDifferencePerturbation);
    } else if (plusValid) {
      column.beDifferenceOf(plusForce, baselineForce);
      column.times(1.0 / finiteDifferencePerturbation);
    } else if (minusValid) {
      column.beDifferenceOf(baselineForce, minusForce);
      column.times(1.0 / finiteDifferencePerturbation);
    } else {
      OOFEM_ERROR(
        "StructuralPenaltyContactBoundaryCondition: both finite-difference "
        "perturbations left the baseline master element");
    }
    numericalTangent.setColumn(column, col);

    dof->updateUnknownsDictionary(tStep, VM_Total, originalUnknown);
    if (!this->updateFiniteDifferencePairGeometry(
          contactPair, masterElementId, featureType, featureIndex, tStep)) {
      OOFEM_ERROR(
        "StructuralPenaltyContactBoundaryCondition: failed to restore baseline "
        "contact projection during finite-difference tangent evaluation");
    }
  }

  if (!this->updateFiniteDifferencePairGeometry(
        contactPair, masterElementId, featureType, featureIndex, tStep)) {
    OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: failed to restore baseline contact projection after finite-difference check");
  }
  finiteDifferenceCheckInProgress = false;
  finiteDifferenceContactPair = nullptr;
  finiteDifferenceContactArea = 0.0;

  return numericalTangent;
}

void
StructuralPenaltyContactBoundaryCondition::dumpFiniteDifferenceKinematicCheck(
    std::ostream &output, const FloatMatrix &suppliedTangent,
    const FloatMatrix &numericalTangent,
    ContactPair *contactPair, TimeStep *tStep)
{
  struct KinematicState {
    FloatArray masterLocalCoordinates;
    FloatArray normal;
    std::vector<FloatArray> tangents;
    FloatMatrix covariantMetric;
    FloatMatrix contravariantMetric;
    FloatMatrix interpolation;
    FloatArray projectionChord;
    double normalGap = 0.0;
    double normalTraction = 0.0;
    double tangentialNorm = 0.0;
    double frictionBound = 0.0;
    double yieldValue = 0.0;
    double contactArea = 0.0;
    double normalPenalty = 0.0;
    double tangentialPenalty = 0.0;
    ContactProcess mode = ContactProcess::None;
    FloatArray transportedPreviousPhysicalTraction;
    FloatArray transportedCovariantTraction;
    FloatArray penaltyIncrementTraction;
    FloatArray trialTraction;
    FloatArray physicalTrialTraction;
    FloatArray returnedTraction;
    FloatArray physicalReturnedTraction;
    FloatArray totalTraction;
    FloatArray residual;
  };

  const int spatialDimension = surface_dimension + 1;
  const std::vector<Node *> residualNodes = contactPair->giveResidualNodes();
  const std::vector<Node *> linearizationNodes = contactPair->giveLinearizationNodes();
  const int nColumns = static_cast<int>(linearizationNodes.size()) * dofs.giveSize();
  if (contactPair->hasMasterFacetTransition()
      || residualNodes != linearizationNodes
      || numericalTangent.giveNumberOfColumns() != nColumns) {
    output << "kinematic_fd skipped_nonlocal_history\n";
    return;
  }

  auto captureState = [&]() {
    KinematicState state;
    auto *masterPoint = dynamic_cast<FEContactPoint *>(contactPair->giveMasterContactPoint());
    if (masterPoint == nullptr) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: kinematic FD audit requires an FE master point");
    }
    state.masterLocalCoordinates = masterPoint->giveLocalCoordinates();
    state.normal = contactPair->giveNormalVector();
    const double normalNorm = state.normal.computeNorm();
    if (normalNorm <= 0.0) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: zero normal in kinematic FD audit");
    }
    state.normal /= normalNorm;
    state.tangents = contactPair->giveTangentVectors();
    state.covariantMetric = this->computeCovariantMetric(state.tangents);
    state.contravariantMetric = this->computeContravariantMetric(state.tangents);
    contactPair->computeNmatrix(state.interpolation);
    state.projectionChord = contactPair->computeContactPointDisplacement(tStep);
    const ContactTractionState tractionState = this->computeTractionState(contactPair, tStep);
    state.normalGap = contactPair->giveNormalGap();
    state.normalTraction = tractionState.normalTraction;
    state.tangentialNorm = tractionState.tangentialNorm;
    state.frictionBound = tractionState.frictionBound;
    state.yieldValue = state.tangentialNorm - state.frictionBound;
    state.contactArea = this->giveContactAreaForLinearization(contactPair, tStep);
    state.normalPenalty = this->giveNormalPenalty(contactPair);
    state.tangentialPenalty = this->giveTangentialPenalty(contactPair);
    state.mode = tractionState.mode;
    state.trialTraction = tractionState.trialTangentialTraction;
    state.returnedTraction = tractionState.tangentialTraction;

    state.transportedPreviousPhysicalTraction.resize(spatialDimension);
    state.transportedPreviousPhysicalTraction.zero();
    if (frictionShouldBeConsidered(friction, contactPair)) {
      state.transportedPreviousPhysicalTraction =
        this->computeTransportedPreviousPhysicalTraction(contactPair);
    }
    state.transportedCovariantTraction.resize(surface_dimension);
    state.penaltyIncrementTraction.resize(surface_dimension);
    for (int i = 0; i < surface_dimension; ++i) {
      state.transportedCovariantTraction.at(i + 1) =
        state.transportedPreviousPhysicalTraction.dotProduct(state.tangents[i]);
      state.penaltyIncrementTraction.at(i + 1) = state.trialTraction.at(i + 1)
        - state.transportedCovariantTraction.at(i + 1);
    }

    state.physicalTrialTraction.resize(spatialDimension);
    state.physicalTrialTraction.zero();
    state.physicalReturnedTraction.resize(spatialDimension);
    state.physicalReturnedTraction.zero();
    for (int i = 0; i < surface_dimension; ++i) {
      for (int j = 0; j < surface_dimension; ++j) {
        state.physicalTrialTraction.add(
          state.contravariantMetric.at(i + 1, j + 1)
            * state.trialTraction.at(j + 1),
          state.tangents[i]);
        state.physicalReturnedTraction.add(
          state.contravariantMetric.at(i + 1, j + 1)
            * state.returnedTraction.at(j + 1),
          state.tangents[i]);
      }
    }
    state.totalTraction = state.normalTraction * state.normal;
    if (frictionShouldBeConsidered(friction, contactPair)) {
      state.totalTraction += state.physicalReturnedTraction;
    }
    state.residual.beTProductOf(state.interpolation, state.totalTraction);
    state.residual.times(state.contactArea);
    return state;
  };

  auto *baselineMasterPoint = dynamic_cast<FEContactPoint *>(contactPair->giveMasterContactPoint());
  if (baselineMasterPoint == nullptr) {
    OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: kinematic FD audit requires an FE master point");
  }
  const int masterElementId = baselineMasterPoint->giveContactElementId();
  const ContactFeatureType featureType = baselineMasterPoint->giveContactFeatureType();
  const int featureIndex = baselineMasterPoint->giveContactFeatureIndex();

  std::vector<Dof *> perturbationDofs;
  std::vector<double> originalUnknowns;
  perturbationDofs.reserve(nColumns);
  originalUnknowns.reserve(nColumns);
  for (Node *node : linearizationNodes) {
    for (int component = 1; component <= dofs.giveSize(); ++component) {
      Dof *dof = node->giveDofWithID(dofs.at(component));
      if (dof == nullptr || !dof->isPrimaryDof()) {
        OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: kinematic FD audit requires primary displacement DOFs");
      }
      perturbationDofs.push_back(dof);
      originalUnknowns.push_back(dof->giveUnknown(VM_Total, tStep));
    }
  }

  FloatMatrix numericalLocal(surface_dimension, nColumns);
  FloatMatrix numericalNormal(spatialDimension, nColumns);
  std::vector<FloatMatrix> numericalTangents(
    surface_dimension, FloatMatrix(spatialDimension, nColumns));
  FloatMatrix numericalCovariant(surface_dimension * surface_dimension, nColumns);
  FloatMatrix numericalContravariant(surface_dimension * surface_dimension, nColumns);
  const int nResidualRows = numericalTangent.giveNumberOfRows();
  FloatMatrix numericalInterpolation(spatialDimension * nResidualRows, nColumns);
  FloatMatrix numericalChord(spatialDimension, nColumns);
  FloatMatrix numericalGap(1, nColumns);
  FloatMatrix numericalNormalTraction(1, nColumns);
  FloatMatrix numericalTangentialNorm(1, nColumns);
  FloatMatrix numericalFrictionBound(1, nColumns);
  FloatMatrix numericalYieldValue(1, nColumns);
  FloatMatrix numericalContactArea(1, nColumns);
  FloatMatrix numericalNormalPenalty(1, nColumns);
  FloatMatrix numericalTangentialPenalty(1, nColumns);
  FloatMatrix numericalTransportedPreviousPhysical(spatialDimension, nColumns);
  FloatMatrix numericalTransportedCovariant(surface_dimension, nColumns);
  FloatMatrix numericalPenaltyIncrement(surface_dimension, nColumns);
  FloatMatrix numericalTrial(surface_dimension, nColumns);
  FloatMatrix numericalPhysicalTrial(spatialDimension, nColumns);
  FloatMatrix numericalReturned(surface_dimension, nColumns);
  FloatMatrix numericalPhysicalReturned(spatialDimension, nColumns);
  FloatMatrix numericalTotalTraction(spatialDimension, nColumns);
  FloatMatrix numericalResidual(nResidualRows, nColumns);
  numericalLocal.zero();
  numericalNormal.zero();
  numericalCovariant.zero();
  numericalContravariant.zero();
  numericalInterpolation.zero();
  numericalChord.zero();
  numericalGap.zero();
  numericalNormalTraction.zero();
  numericalTangentialNorm.zero();
  numericalFrictionBound.zero();
  numericalYieldValue.zero();
  numericalContactArea.zero();
  numericalNormalPenalty.zero();
  numericalTangentialPenalty.zero();
  numericalTransportedPreviousPhysical.zero();
  numericalTransportedCovariant.zero();
  numericalPenaltyIncrement.zero();
  numericalTrial.zero();
  numericalPhysicalTrial.zero();
  numericalReturned.zero();
  numericalPhysicalReturned.zero();
  numericalTotalTraction.zero();
  numericalResidual.zero();
  for (FloatMatrix &matrix : numericalTangents) {
    matrix.zero();
  }

  finiteDifferenceContactPair = contactPair;
  finiteDifferenceContactArea = contactPair->giveContactArea(tStep);
  finiteDifferenceCheckInProgress = true;
  const KinematicState baseline = captureState();
  bool returnMapBranchChanged = false;

  for (int col = 1; col <= nColumns; ++col) {
    Dof *dof = perturbationDofs[col - 1];
    const double originalUnknown = originalUnknowns[col - 1];
    dof->updateUnknownsDictionary(
      tStep, VM_Total, originalUnknown + finiteDifferencePerturbation);
    const bool plusValid = this->updateFiniteDifferencePairGeometry(
      contactPair, masterElementId, featureType, featureIndex, tStep);
    KinematicState plus;
    if (plusValid) {
      plus = captureState();
    }

    dof->updateUnknownsDictionary(
      tStep, VM_Total, originalUnknown - finiteDifferencePerturbation);
    const bool minusValid = this->updateFiniteDifferencePairGeometry(
      contactPair, masterElementId, featureType, featureIndex, tStep);
    KinematicState minus;
    if (minusValid) {
      minus = captureState();
    }
    if (!plusValid || !minusValid) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: kinematic FD audit left the baseline master element");
    }
    returnMapBranchChanged = returnMapBranchChanged
      || plus.mode != baseline.mode || minus.mode != baseline.mode;

    const double scale = 0.5 / finiteDifferencePerturbation;
    numericalGap.at(1, col) = scale * (plus.normalGap - minus.normalGap);
    numericalNormalTraction.at(1, col) = scale
      * (plus.normalTraction - minus.normalTraction);
    numericalTangentialNorm.at(1, col) = scale
      * (plus.tangentialNorm - minus.tangentialNorm);
    numericalFrictionBound.at(1, col) = scale
      * (plus.frictionBound - minus.frictionBound);
    numericalYieldValue.at(1, col) = scale
      * (plus.yieldValue - minus.yieldValue);
    numericalContactArea.at(1, col) = scale
      * (plus.contactArea - minus.contactArea);
    numericalNormalPenalty.at(1, col) = scale
      * (plus.normalPenalty - minus.normalPenalty);
    numericalTangentialPenalty.at(1, col) = scale
      * (plus.tangentialPenalty - minus.tangentialPenalty);
    for (int i = 1; i <= surface_dimension; ++i) {
      numericalLocal.at(i, col) = scale
        * (plus.masterLocalCoordinates.at(i) - minus.masterLocalCoordinates.at(i));
      numericalTransportedCovariant.at(i, col) = scale
        * (plus.transportedCovariantTraction.at(i)
           - minus.transportedCovariantTraction.at(i));
      numericalPenaltyIncrement.at(i, col) = scale
        * (plus.penaltyIncrementTraction.at(i)
           - minus.penaltyIncrementTraction.at(i));
      numericalTrial.at(i, col) = scale
        * (plus.trialTraction.at(i) - minus.trialTraction.at(i));
      numericalReturned.at(i, col) = scale
        * (plus.returnedTraction.at(i) - minus.returnedTraction.at(i));
    }
    for (int r = 1; r <= spatialDimension; ++r) {
      numericalNormal.at(r, col) = scale
        * (plus.normal.at(r) - minus.normal.at(r));
      numericalChord.at(r, col) = scale
        * (plus.projectionChord.at(r) - minus.projectionChord.at(r));
      numericalPhysicalTrial.at(r, col) = scale
        * (plus.physicalTrialTraction.at(r) - minus.physicalTrialTraction.at(r));
      numericalTransportedPreviousPhysical.at(r, col) = scale
        * (plus.transportedPreviousPhysicalTraction.at(r)
           - minus.transportedPreviousPhysicalTraction.at(r));
      numericalPhysicalReturned.at(r, col) = scale
        * (plus.physicalReturnedTraction.at(r)
           - minus.physicalReturnedTraction.at(r));
      numericalTotalTraction.at(r, col) = scale
        * (plus.totalTraction.at(r) - minus.totalTraction.at(r));
      for (int i = 0; i < surface_dimension; ++i) {
        numericalTangents[i].at(r, col) = scale
          * (plus.tangents[i].at(r) - minus.tangents[i].at(r));
      }
      for (int residualRow = 1; residualRow <= nResidualRows; ++residualRow) {
        const int interpolationRow = (residualRow - 1) * spatialDimension + r;
        numericalInterpolation.at(interpolationRow, col) = scale
          * (plus.interpolation.at(r, residualRow)
             - minus.interpolation.at(r, residualRow));
      }
    }
    for (int row = 1; row <= nResidualRows; ++row) {
      numericalResidual.at(row, col) = scale
        * (plus.residual.at(row) - minus.residual.at(row));
    }
    for (int i = 1; i <= surface_dimension; ++i) {
      for (int j = 1; j <= surface_dimension; ++j) {
        const int row = (i - 1) * surface_dimension + j;
        numericalCovariant.at(row, col) = scale
          * (plus.covariantMetric.at(i, j) - minus.covariantMetric.at(i, j));
        numericalContravariant.at(row, col) = scale
          * (plus.contravariantMetric.at(i, j) - minus.contravariantMetric.at(i, j));
      }
    }

    dof->updateUnknownsDictionary(tStep, VM_Total, originalUnknown);
    if (!this->updateFiniteDifferencePairGeometry(
          contactPair, masterElementId, featureType, featureIndex, tStep)) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: failed to restore geometry in kinematic FD audit");
    }
  }

  finiteDifferenceCheckInProgress = false;
  finiteDifferenceContactPair = nullptr;
  finiteDifferenceContactArea = 0.0;

  FloatMatrix N;
  contactPair->computeNmatrix(N);
  std::vector<FloatMatrix> dNs;
  contactPair->compute_dNdxi_matrices(dNs);
  if (N.giveNumberOfColumns() != nColumns
      || static_cast<int>(dNs.size()) != surface_dimension) {
    output << "kinematic_fd skipped_column_layout\n";
    return;
  }

  FloatMatrix curvature;
  contactPair->computeCurvature(curvature, tStep);
  std::vector<std::vector<FloatArray>> secondBaseVectors;
  contactPair->computeSecondBaseVectors(secondBaseVectors, tStep);
  FloatMatrix projectionHessian(surface_dimension, surface_dimension);
  for (int i = 1; i <= surface_dimension; ++i) {
    for (int j = 1; j <= surface_dimension; ++j) {
      projectionHessian.at(i, j) = baseline.covariantMetric.at(i, j)
        - contactPair->giveNormalGap() * curvature.at(i, j);
    }
  }
  FloatMatrix inverseProjectionHessian;
  inverseProjectionHessian.beInverseOf(projectionHessian);

  FloatMatrix localRightHandSide(surface_dimension, nColumns);
  for (int i = 0; i < surface_dimension; ++i) {
    for (int col = 1; col <= nColumns; ++col) {
      double value = 0.0;
      double normalDerivative = 0.0;
      for (int r = 1; r <= spatialDimension; ++r) {
        value += baseline.tangents[i].at(r) * N.at(r, col);
        normalDerivative += baseline.normal.at(r) * dNs[i].at(r, col);
      }
      localRightHandSide.at(i + 1, col) = value
        - contactPair->giveNormalGap() * normalDerivative;
    }
  }
  FloatMatrix analyticalLocal;
  analyticalLocal.beProductOf(inverseProjectionHessian, localRightHandSide);

  std::vector<FloatMatrix> analyticalTangents(
    surface_dimension, FloatMatrix(spatialDimension, nColumns));
  for (int i = 0; i < surface_dimension; ++i) {
    for (int r = 1; r <= spatialDimension; ++r) {
      for (int col = 1; col <= nColumns; ++col) {
        double value = -dNs[i].at(r, col);
        for (int j = 0; j < surface_dimension; ++j) {
          value += secondBaseVectors[i][j].at(r)
            * analyticalLocal.at(j + 1, col);
        }
        analyticalTangents[i].at(r, col) = value;
      }
    }
  }

  std::vector<FloatArray> contravariantTangents(surface_dimension);
  for (int i = 0; i < surface_dimension; ++i) {
    contravariantTangents[i].resize(spatialDimension);
    contravariantTangents[i].zero();
    for (int j = 0; j < surface_dimension; ++j) {
      contravariantTangents[i].add(
        baseline.contravariantMetric.at(i + 1, j + 1), baseline.tangents[j]);
    }
  }
  FloatMatrix analyticalNormal(spatialDimension, nColumns);
  analyticalNormal.zero();
  for (int col = 1; col <= nColumns; ++col) {
    for (int i = 0; i < surface_dimension; ++i) {
      double normalDotTangentDerivative = 0.0;
      for (int r = 1; r <= spatialDimension; ++r) {
        normalDotTangentDerivative += baseline.normal.at(r)
          * analyticalTangents[i].at(r, col);
      }
      for (int r = 1; r <= spatialDimension; ++r) {
        analyticalNormal.at(r, col) -= contravariantTangents[i].at(r)
          * normalDotTangentDerivative;
      }
    }
  }

  FloatMatrix analyticalCovariant(surface_dimension * surface_dimension, nColumns);
  FloatMatrix analyticalContravariant(surface_dimension * surface_dimension, nColumns);
  for (int i = 0; i < surface_dimension; ++i) {
    for (int j = 0; j < surface_dimension; ++j) {
      const int row = i * surface_dimension + j + 1;
      for (int col = 1; col <= nColumns; ++col) {
        double value = 0.0;
        for (int r = 1; r <= spatialDimension; ++r) {
          value += analyticalTangents[i].at(r, col) * baseline.tangents[j].at(r)
            + baseline.tangents[i].at(r) * analyticalTangents[j].at(r, col);
        }
        analyticalCovariant.at(row, col) = value;
      }
    }
  }
  for (int i = 0; i < surface_dimension; ++i) {
    for (int j = 0; j < surface_dimension; ++j) {
      const int row = i * surface_dimension + j + 1;
      for (int col = 1; col <= nColumns; ++col) {
        double value = 0.0;
        for (int k = 0; k < surface_dimension; ++k) {
          for (int l = 0; l < surface_dimension; ++l) {
            value -= baseline.contravariantMetric.at(i + 1, k + 1)
              * analyticalCovariant.at(k * surface_dimension + l + 1, col)
              * baseline.contravariantMetric.at(l + 1, j + 1);
          }
        }
        analyticalContravariant.at(row, col) = value;
      }
    }
  }

  // Updating a pair replaces its FEContactPoint, so the pointer captured before
  // the perturbation loop is no longer valid even though the baseline geometry
  // has been restored.  Reacquire it before using element/interpolation data.
  auto *restoredMasterPoint =
    dynamic_cast<FEContactPoint *>(contactPair->giveMasterContactPoint());
  if (restoredMasterPoint == nullptr) {
    OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: failed to restore FE master point after kinematic FD audit");
  }
  FloatMatrix currentMasterN;
  restoredMasterPoint->computeNmatrix(currentMasterN);
  ContactElement *previousMaster = contactPair->givePreviousMasterContactElement();
  if (previousMaster == nullptr
      || previousMaster != restoredMasterPoint->giveContactElement()) {
    output << "kinematic_fd skipped_previous_master\n";
    return;
  }
  FloatMatrix previousMasterN;
  previousMaster->computeNmatrixAt(
    contactPair->givePreviousMasterLocalCoordinates(), previousMasterN);
  if (currentMasterN.giveNumberOfColumns() != previousMasterN.giveNumberOfColumns()) {
    output << "kinematic_fd skipped_previous_layout\n";
    return;
  }
  FloatMatrix analyticalChord(spatialDimension, nColumns);
  analyticalChord.zero();
  const int masterColumns = currentMasterN.giveNumberOfColumns();
  for (int r = 1; r <= spatialDimension; ++r) {
    for (int col = 1; col <= masterColumns; ++col) {
      analyticalChord.at(r, col) = currentMasterN.at(r, col)
        - previousMasterN.at(r, col);
    }
    for (int col = 1; col <= nColumns; ++col) {
      for (int i = 0; i < surface_dimension; ++i) {
        analyticalChord.at(r, col) += baseline.tangents[i].at(r)
          * analyticalLocal.at(i + 1, col);
      }
    }
  }

  FloatMatrix analyticalTransportedCovariant(surface_dimension, nColumns);
  FloatMatrix analyticalPenaltyIncrement(surface_dimension, nColumns);
  FloatMatrix analyticalTrial(surface_dimension, nColumns);
  analyticalTransportedCovariant.zero();
  analyticalPenaltyIncrement.zero();
  analyticalTrial.zero();
  const double tangentialPenalty = this->giveTangentialPenalty(contactPair);
  for (int i = 0; i < surface_dimension; ++i) {
    for (int col = 1; col <= nColumns; ++col) {
      double tangentDerivative = 0.0;
      double chordDerivative = 0.0;
      for (int r = 1; r <= spatialDimension; ++r) {
        tangentDerivative += baseline.transportedPreviousPhysicalTraction.at(r)
          * analyticalTangents[i].at(r, col);
        chordDerivative += analyticalChord.at(r, col) * baseline.tangents[i].at(r)
          + baseline.projectionChord.at(r) * analyticalTangents[i].at(r, col);
      }
      analyticalTransportedCovariant.at(i + 1, col) = tangentDerivative;
      analyticalPenaltyIncrement.at(i + 1, col) =
        -tangentialPenalty * chordDerivative;
      analyticalTrial.at(i + 1, col) =
        analyticalTransportedCovariant.at(i + 1, col)
        + analyticalPenaltyIncrement.at(i + 1, col);
    }
  }

  FloatMatrix analyticalPhysicalTrial(spatialDimension, nColumns);
  analyticalPhysicalTrial.zero();
  for (int r = 1; r <= spatialDimension; ++r) {
    for (int col = 1; col <= nColumns; ++col) {
      double value = 0.0;
      for (int i = 0; i < surface_dimension; ++i) {
        for (int j = 0; j < surface_dimension; ++j) {
          const double aij = baseline.contravariantMetric.at(i + 1, j + 1);
          const double trialJ = baseline.trialTraction.at(j + 1);
          value += analyticalTangents[i].at(r, col) * aij * trialJ;
          value += baseline.tangents[i].at(r)
            * analyticalContravariant.at(i * surface_dimension + j + 1, col)
            * trialJ;
          value += baseline.tangents[i].at(r) * aij
            * analyticalTrial.at(j + 1, col);
        }
      }
      analyticalPhysicalTrial.at(r, col) = value;
    }
  }

  FloatMatrix analyticalInterpolation(spatialDimension * nResidualRows, nColumns);
  analyticalInterpolation.zero();
  for (int residualRow = 1; residualRow <= nResidualRows; ++residualRow) {
    for (int r = 1; r <= spatialDimension; ++r) {
      const int interpolationRow = (residualRow - 1) * spatialDimension + r;
      for (int col = 1; col <= nColumns; ++col) {
        for (int i = 0; i < surface_dimension; ++i) {
          analyticalInterpolation.at(interpolationRow, col) +=
            dNs[i].at(r, residualRow) * analyticalLocal.at(i + 1, col);
        }
      }
    }
  }

  FloatMatrix analyticalGap(1, nColumns);
  FloatMatrix analyticalNormalTraction(1, nColumns);
  FloatMatrix analyticalTangentialNorm(1, nColumns);
  FloatMatrix analyticalFrictionBound(1, nColumns);
  FloatMatrix analyticalYieldValue(1, nColumns);
  FloatMatrix analyticalContactArea(1, nColumns);
  FloatMatrix analyticalNormalPenalty(1, nColumns);
  FloatMatrix analyticalTangentialPenalty(1, nColumns);
  analyticalTangentialNorm.zero();
  analyticalFrictionBound.zero();
  analyticalYieldValue.zero();
  analyticalContactArea.zero();
  analyticalNormalPenalty.zero();
  analyticalTangentialPenalty.zero();
  const double normalPenalty = this->giveNormalPenalty(contactPair);
  for (int col = 1; col <= nColumns; ++col) {
    double gapDerivative = 0.0;
    for (int r = 1; r <= spatialDimension; ++r) {
      gapDerivative += baseline.normal.at(r) * N.at(r, col);
    }
    analyticalGap.at(1, col) = gapDerivative;
    analyticalNormalTraction.at(1, col) = normalPenalty * gapDerivative;
  }

  FloatMatrix analyticalTransportedPreviousPhysical(spatialDimension, nColumns);
  analyticalTransportedPreviousPhysical.zero();
  FloatMatrix analyticalReturned(surface_dimension, nColumns);
  analyticalReturned.zero();
  if (frictionShouldBeConsidered(friction, contactPair)) {
    if (baseline.tangentialNorm > 0.0) {
      for (int col = 1; col <= nColumns; ++col) {
        double squaredNormDerivative = 0.0;
        for (int i = 0; i < surface_dimension; ++i) {
          for (int j = 0; j < surface_dimension; ++j) {
            const double aij = baseline.contravariantMetric.at(i + 1, j + 1);
            const double trialI = baseline.trialTraction.at(i + 1);
            const double trialJ = baseline.trialTraction.at(j + 1);
            squaredNormDerivative += analyticalTrial.at(i + 1, col) * aij * trialJ;
            squaredNormDerivative += trialI
              * analyticalContravariant.at(i * surface_dimension + j + 1, col)
              * trialJ;
            squaredNormDerivative += trialI * aij * analyticalTrial.at(j + 1, col);
          }
        }
        analyticalTangentialNorm.at(1, col) =
          0.5 * squaredNormDerivative / baseline.tangentialNorm;
      }
    }
    const double normalTractionSign =
      baseline.normalTraction < 0.0 ? -1.0 : (baseline.normalTraction > 0.0 ? 1.0 : 0.0);
    for (int col = 1; col <= nColumns; ++col) {
      analyticalFrictionBound.at(1, col) = friction * normalTractionSign
        * analyticalNormalTraction.at(1, col);
      analyticalYieldValue.at(1, col) = analyticalTangentialNorm.at(1, col)
        - analyticalFrictionBound.at(1, col);
    }

    if (baseline.tangentialNorm > 0.0) {
      double returnedMagnitude = 0.0;
      double derivativeWrtTrial = 0.0;
      double derivativeWrtBound = 0.0;
      this->computeFrictionReturnMagnitude(
        returnedMagnitude, derivativeWrtTrial, derivativeWrtBound,
        baseline.tangentialNorm, baseline.frictionBound);
      const double returnScale = returnedMagnitude / baseline.tangentialNorm;
      for (int i = 1; i <= surface_dimension; ++i) {
        for (int col = 1; col <= nColumns; ++col) {
          const double scaleDerivative =
            (derivativeWrtTrial * baseline.tangentialNorm - returnedMagnitude)
              / (baseline.tangentialNorm * baseline.tangentialNorm)
              * analyticalTangentialNorm.at(1, col)
            + derivativeWrtBound / baseline.tangentialNorm
              * analyticalFrictionBound.at(1, col);
          analyticalReturned.at(i, col) = returnScale * analyticalTrial.at(i, col)
            + baseline.trialTraction.at(i) * scaleDerivative;
        }
      }
    }
  }

  FloatMatrix analyticalPhysicalReturned(spatialDimension, nColumns);
  analyticalPhysicalReturned.zero();
  for (int r = 1; r <= spatialDimension; ++r) {
    for (int col = 1; col <= nColumns; ++col) {
      for (int i = 0; i < surface_dimension; ++i) {
        for (int j = 0; j < surface_dimension; ++j) {
          const double aij = baseline.contravariantMetric.at(i + 1, j + 1);
          const double returnedJ = baseline.returnedTraction.at(j + 1);
          analyticalPhysicalReturned.at(r, col) +=
            analyticalTangents[i].at(r, col) * aij * returnedJ;
          analyticalPhysicalReturned.at(r, col) += baseline.tangents[i].at(r)
            * analyticalContravariant.at(i * surface_dimension + j + 1, col)
            * returnedJ;
          analyticalPhysicalReturned.at(r, col) += baseline.tangents[i].at(r)
            * aij * analyticalReturned.at(j + 1, col);
        }
      }
    }
  }

  FloatMatrix analyticalTotalTraction(spatialDimension, nColumns);
  analyticalTotalTraction.zero();
  for (int r = 1; r <= spatialDimension; ++r) {
    for (int col = 1; col <= nColumns; ++col) {
      analyticalTotalTraction.at(r, col) =
        analyticalNormalTraction.at(1, col) * baseline.normal.at(r)
        + baseline.normalTraction * analyticalNormal.at(r, col);
      if (frictionShouldBeConsidered(friction, contactPair)) {
        analyticalTotalTraction.at(r, col) +=
          analyticalPhysicalReturned.at(r, col);
      }
    }
  }

  FloatMatrix analyticalResidualInterpolation(nResidualRows, nColumns);
  FloatMatrix numericalResidualInterpolation(nResidualRows, nColumns);
  FloatMatrix analyticalResidualTraction(nResidualRows, nColumns);
  FloatMatrix numericalResidualTraction(nResidualRows, nColumns);
  FloatMatrix analyticalResidualArea(nResidualRows, nColumns);
  FloatMatrix numericalResidualArea(nResidualRows, nColumns);
  analyticalResidualInterpolation.zero();
  numericalResidualInterpolation.zero();
  analyticalResidualTraction.zero();
  numericalResidualTraction.zero();
  analyticalResidualArea.zero();
  numericalResidualArea.zero();
  for (int row = 1; row <= nResidualRows; ++row) {
    for (int col = 1; col <= nColumns; ++col) {
      double baselineUnscaledResidual = 0.0;
      for (int r = 1; r <= spatialDimension; ++r) {
        const int interpolationRow = (row - 1) * spatialDimension + r;
        analyticalResidualInterpolation.at(row, col) +=
          analyticalInterpolation.at(interpolationRow, col)
          * baseline.totalTraction.at(r);
        numericalResidualInterpolation.at(row, col) +=
          numericalInterpolation.at(interpolationRow, col)
          * baseline.totalTraction.at(r);
        analyticalResidualTraction.at(row, col) +=
          N.at(r, row) * analyticalTotalTraction.at(r, col);
        numericalResidualTraction.at(row, col) +=
          N.at(r, row) * numericalTotalTraction.at(r, col);
        baselineUnscaledResidual += N.at(r, row) * baseline.totalTraction.at(r);
      }
      analyticalResidualArea.at(row, col) = baselineUnscaledResidual
        * analyticalContactArea.at(1, col);
      numericalResidualArea.at(row, col) = baselineUnscaledResidual
        * numericalContactArea.at(1, col);
    }
  }
  analyticalResidualInterpolation.times(baseline.contactArea);
  numericalResidualInterpolation.times(baseline.contactArea);
  analyticalResidualTraction.times(baseline.contactArea);
  numericalResidualTraction.times(baseline.contactArea);
  FloatMatrix analyticalResidual = analyticalResidualInterpolation;
  analyticalResidual += analyticalResidualTraction;
  analyticalResidual += analyticalResidualArea;

  auto reportError = [&](const char *label, const FloatMatrix &analytical,
                         const FloatMatrix &numerical) {
    double analyticalNormSquared = 0.0;
    double differenceNormSquared = 0.0;
    for (int i = 1; i <= analytical.giveNumberOfRows(); ++i) {
      for (int j = 1; j <= analytical.giveNumberOfColumns(); ++j) {
        analyticalNormSquared += analytical.at(i, j) * analytical.at(i, j);
        const double difference = analytical.at(i, j) - numerical.at(i, j);
        differenceNormSquared += difference * difference;
      }
    }
    const double analyticalNorm = std::sqrt(analyticalNormSquared);
    const double differenceNorm = std::sqrt(differenceNormSquared);
    output << "kinematic_fd " << label
           << " relative_error " << differenceNorm / std::max(analyticalNorm, 1.e-30)
           << " analytical_norm " << analyticalNorm
           << " difference_norm " << differenceNorm << '\n';
  };

  writeVectorBlock(output, "kinematic_baseline_master_local", baseline.masterLocalCoordinates);
  writeVectorBlock(output, "kinematic_baseline_normal", baseline.normal);
  for (int i = 0; i < surface_dimension; ++i) {
    const std::string label = "kinematic_baseline_tangent_" + std::to_string(i + 1);
    writeVectorBlock(output, label.c_str(), baseline.tangents[i]);
  }
  writeMatrixBlock(output, "kinematic_baseline_covariant_metric", baseline.covariantMetric);
  writeMatrixBlock(output, "kinematic_baseline_contravariant_metric", baseline.contravariantMetric);
  writeMatrixBlock(output, "kinematic_baseline_interpolation", baseline.interpolation);
  writeVectorBlock(output, "kinematic_baseline_projection_chord", baseline.projectionChord);
  output << "kinematic_baseline_normal_gap " << baseline.normalGap << '\n';
  output << "kinematic_baseline_normal_traction " << baseline.normalTraction << '\n';
  output << "kinematic_baseline_tangential_norm " << baseline.tangentialNorm << '\n';
  output << "kinematic_baseline_friction_bound " << baseline.frictionBound << '\n';
  output << "kinematic_baseline_yield_value " << baseline.yieldValue << '\n';
  output << "kinematic_baseline_contact_area " << baseline.contactArea << '\n';
  output << "kinematic_baseline_normal_penalty " << baseline.normalPenalty << '\n';
  output << "kinematic_baseline_tangential_penalty " << baseline.tangentialPenalty << '\n';
  writeVectorBlock(output, "kinematic_baseline_transported_previous_physical_traction",
                   baseline.transportedPreviousPhysicalTraction);
  writeVectorBlock(output, "kinematic_baseline_transported_covariant_traction",
                   baseline.transportedCovariantTraction);
  writeVectorBlock(output, "kinematic_baseline_penalty_increment_traction",
                   baseline.penaltyIncrementTraction);
  writeVectorBlock(output, "kinematic_baseline_trial_traction", baseline.trialTraction);
  writeVectorBlock(output, "kinematic_baseline_physical_trial_traction", baseline.physicalTrialTraction);
  writeVectorBlock(output, "kinematic_baseline_returned_traction", baseline.returnedTraction);
  writeVectorBlock(output, "kinematic_baseline_physical_returned_traction",
                   baseline.physicalReturnedTraction);
  writeVectorBlock(output, "kinematic_baseline_total_traction", baseline.totalTraction);
  writeVectorBlock(output, "kinematic_baseline_residual", baseline.residual);
  reportError("master_local", analyticalLocal, numericalLocal);
  for (int i = 0; i < surface_dimension; ++i) {
    const std::string label = "tangent_" + std::to_string(i + 1);
    reportError(label.c_str(), analyticalTangents[i], numericalTangents[i]);
  }
  reportError("normal", analyticalNormal, numericalNormal);
  reportError("covariant_metric", analyticalCovariant, numericalCovariant);
  reportError("contravariant_metric", analyticalContravariant, numericalContravariant);
  reportError("interpolation", analyticalInterpolation, numericalInterpolation);
  reportError("projection_chord", analyticalChord, numericalChord);
  reportError("normal_gap", analyticalGap, numericalGap);
  reportError("normal_traction", analyticalNormalTraction, numericalNormalTraction);
  reportError("contact_area", analyticalContactArea, numericalContactArea);
  reportError("normal_penalty", analyticalNormalPenalty, numericalNormalPenalty);
  reportError("tangential_penalty", analyticalTangentialPenalty, numericalTangentialPenalty);

  if (frictionShouldBeConsidered(friction, contactPair)) {
    reportError("transported_previous_physical_traction",
                analyticalTransportedPreviousPhysical,
                numericalTransportedPreviousPhysical);
    reportError("transported_covariant_traction",
                analyticalTransportedCovariant, numericalTransportedCovariant);
    reportError("penalty_increment_traction",
                analyticalPenaltyIncrement, numericalPenaltyIncrement);
    reportError("trial_traction", analyticalTrial, numericalTrial);
    reportError("physical_trial_traction", analyticalPhysicalTrial, numericalPhysicalTrial);
    reportError("tangential_norm", analyticalTangentialNorm, numericalTangentialNorm);
    reportError("friction_bound", analyticalFrictionBound, numericalFrictionBound);
    reportError("yield_value", analyticalYieldValue, numericalYieldValue);
    if (returnMapBranchChanged) {
      output << "kinematic_fd return_map_branch_changed\n";
    } else {
      reportError("returned_traction", analyticalReturned, numericalReturned);
      reportError("physical_returned_traction",
                  analyticalPhysicalReturned, numericalPhysicalReturned);
    }
  }

  if (!returnMapBranchChanged) {
    reportError("total_traction", analyticalTotalTraction, numericalTotalTraction);
    reportError("residual_interpolation_contribution",
                analyticalResidualInterpolation, numericalResidualInterpolation);
    reportError("residual_traction_contribution",
                analyticalResidualTraction, numericalResidualTraction);
    reportError("residual_area_contribution",
                analyticalResidualArea, numericalResidualArea);
    reportError("residual_capture", analyticalResidual, numericalResidual);
    reportError("residual_capture_vs_contact_fd", numericalResidual, numericalTangent);
    reportError("exact_discrete_residual", analyticalResidual, numericalTangent);
  }
  reportError("supplied_residual", numericalTangent, suppliedTangent);
}

void
StructuralPenaltyContactBoundaryCondition::dumpFiniteDifferenceCheck(const FloatMatrix &analyticalTangent, ContactPair *contactPair, TimeStep *tStep)
{
  const std::size_t pairIndex = this->findContactPairIndex(contactPair);
  const FloatMatrix numericalTangent = this->computeFiniteDifferenceTangent(pairIndex, tStep);

  FloatMatrix difference = analyticalTangent;
  difference.subtract(numericalTangent);

  FloatArray residual;
  this->computeInternalForcesFromContact(residual, getContactPairs()[pairIndex].get(), tStep);

  const ContactTractionState tractionState = this->computeTractionState(getContactPairs()[pairIndex].get(), tStep);
  std::ostringstream fileName;
  fileName << finiteDifferenceOutputPrefix
           << "_bc" << this->giveNumber()
           << "_step" << tStep->giveNumber()
           << "_iter" << tStep->giveSubStepNumber()
           << "_pair" << pairIndex
           << ".txt";

  std::ofstream output(fileName.str());
  if (!output.is_open()) {
    OOFEM_WARNING("StructuralPenaltyContactBoundaryCondition: failed to open finite-difference output file '%s'", fileName.str().c_str());
    return;
  }

  output << "contact_fd_check\n";
  output << "bc " << this->giveNumber() << '\n';
  output << "pair " << pairIndex << '\n';
  output << "step " << tStep->giveNumber() << '\n';
  output << "substep " << tStep->giveSubStepNumber() << '\n';
  output << "state_counter " << tStep->giveSolutionStateCounter() << '\n';
  output << "gap " << getContactPairs()[pairIndex]->giveNormalGap() << '\n';
  output << "normal_traction " << tractionState.normalTraction << '\n';
  output << "tangential_norm " << tractionState.tangentialNorm << '\n';
  output << "friction_bound " << tractionState.frictionBound << '\n';
  output << "accumulated_plastic_slip "
         << tractionState.accumulatedPlasticSlip << '\n';
  output << "mode " << static_cast<int>(tractionState.mode) << '\n';
  output << "perturbation " << finiteDifferencePerturbation << '\n';
  output << "normal_reference Konyukhov-Schweizerhof normal tangent\n";
  output << "friction_reference Konyukhov-Schweizerhof return-map friction tangent\n";
  writeVectorBlock(output, "residual", residual);
  writeMatrixBlock(output, "analytical", analyticalTangent);
  writeMatrixBlock(output, "numerical", numericalTangent);
  double analyticalNormSquared = 0.0;
  double differenceNormSquared = 0.0;
  for (int i = 1; i <= analyticalTangent.giveNumberOfRows(); ++i) {
    for (int j = 1; j <= analyticalTangent.giveNumberOfColumns(); ++j) {
      analyticalNormSquared += analyticalTangent.at(i, j) * analyticalTangent.at(i, j);
      differenceNormSquared += difference.at(i, j) * difference.at(i, j);
    }
  }
  const double relativeError = std::sqrt(differenceNormSquared)
    / std::max(std::sqrt(analyticalNormSquared), 1.e-30);
  output << "relative_error " << relativeError << "\n";
  writeMatrixBlock(output, "difference", difference);
  this->dumpFiniteDifferenceKinematicCheck(
    output, analyticalTangent, numericalTangent, contactPair, tStep);

  // gap == 0 and phi == 0 are nonsmooth active-set/return-map boundaries.  A
  // centered difference then averages different one-sided derivatives and
  // cannot equal the selected active branch tangent.  Assert only when both
  // the normal and Coulomb branches are stable under the perturbation.
  const double returnMapPerturbationScale = 10.0 * finiteDifferencePerturbation
    * std::max(this->giveTangentialPenalty(getContactPairs()[pairIndex].get()),
               friction * this->giveNormalPenalty(getContactPairs()[pairIndex].get()));
  const bool stableReturnMapBranch =
    frictionTransition > 0.0
    || std::abs(tractionState.tangentialNorm - tractionState.frictionBound)
         > returnMapPerturbationScale;
  output << "return_map_branch_stable " << stableReturnMapBranch << '\n';
  output.close();
  if (finiteDifferenceRelativeTolerance > 0.0
      && std::abs(getContactPairs()[pairIndex]->giveNormalGap()) > 10.0 * finiteDifferencePerturbation
      && tractionState.mode != ContactProcess::None
      && stableReturnMapBranch
      && relativeError > finiteDifferenceRelativeTolerance) {
    OOFEM_ERROR("Contact tangent finite-difference relative error %e exceeds tolerance %e",
                relativeError, finiteDifferenceRelativeTolerance);
  }
}

ContactTractionState
StructuralPenaltyContactBoundaryCondition :: computeTractionState(ContactPair* contactPair, TimeStep *tStep)
{
  ContactTractionState state;
  state.tangentialTraction.resize(surface_dimension);
  state.tangentialTraction.zero();
  state.trialTangentialTraction.resize(surface_dimension);
  state.trialTangentialTraction.zero();

  const double normal_gap = contactPair->giveNormalGap();
  if (normal_gap > 0.0) {
    return state;
  }

  // Konyukhov-Schweizerhof (2013), Eq. (6.6), p. 136.
  const double tangentialPenalty = this->giveTangentialPenalty(contactPair);
  state.normalTraction = this->giveNormalPenalty(contactPair) * normal_gap;
  if (!frictionShouldBeConsidered(friction, contactPair)) {
    return state;
  }

  // Everything below this point is the experimental, development-phase
  // friction branch (see the @warning on the class declaration) — not yet
  // fully verified. Only reached when `friction` is set nonzero.
  // Objective discrete trial update, Eqs. (6.23)-(6.25), pp. 141-142.
  const auto tangentVectors = contactPair->giveTangentVectors();
  FloatMatrix contravariant_metric = this->computeContravariantMetric(tangentVectors);
  FloatMatrix covariant_metric = this->computeCovariantMetric(tangentVectors);
  const FloatArray transportedPreviousPhysical =
    this->computeTransportedPreviousPhysicalTraction(contactPair);

  FloatArray delta_rho = contactPair->computeContactPointDisplacement(tStep);
  for (int i = 0; i < surface_dimension; i++) {
    const auto& rho_i = tangentVectors[i];
    state.tangentialTraction(i) =
      transportedPreviousPhysical.dotProduct(rho_i);
    for (int j = 0; j < surface_dimension; j++) {
      const double a_ij = covariant_metric(i, j);
      double delta_xi_j = 0.0;
      for (int k = 0; k < surface_dimension; k++) {
        const double a_jk = contravariant_metric(j, k);
        const FloatArray& rho_k = tangentVectors[k];
        delta_xi_j += delta_rho.dotProduct(rho_k) * a_jk;
      }
      // computeContactPointDisplacement is the motion of the master
      // projection in the current configuration. With OOFEM's residual
      // convention (master interpolation enters with a minus sign), the
      // elastic friction traction must have the same sign as this chord.
      // The previous minus sign made the physical friction force assist
      // sliding; an inclined-plane balance exposed that reversal directly.
      state.tangentialTraction(i) += tangentialPenalty * delta_xi_j * a_ij;
    }
  }

  state.trialTangentialTraction = state.tangentialTraction;
  state.tangentialNorm = this->computeTangentialNorm(state.tangentialTraction, contravariant_metric);
  state.accumulatedPlasticSlip = contactPair->giveAccumulatedPlasticSlip();
  const double hardeningModulus = frictionHardening > 0.0
    ? tangentialPenalty * frictionHardening / (1.0 - frictionHardening)
    : 0.0;
  state.frictionBound = friction * std::abs(state.normalTraction)
    + hardeningModulus * state.accumulatedPlasticSlip;
  const double yieldValue = state.tangentialNorm - state.frictionBound;
  // Coulomb return map, Eq. (6.26), p. 142.
  state.mode = yieldValue > 0.0
    ? ContactProcess::Sliding : ContactProcess::Sticking;
  if (state.tangentialNorm > 0.0) {
    // Preserve the original sharp, perfectly plastic Coulomb path exactly
    // when both optional regularizations are disabled.  Besides avoiding
    // needless arithmetic in sticking, this keeps the branch decision tied
    // to yieldValue instead of recomputing it through r/b at the kink.
    if (frictionTransition == 0.0 && frictionHardening == 0.0) {
      if (state.mode == ContactProcess::Sliding) {
        state.tangentialTraction *=
          state.frictionBound / state.tangentialNorm;
      }
    } else {
      double returnedMagnitude = 0.0;
      double unusedTrialDerivative = 0.0;
      double unusedBoundDerivative = 0.0;
      this->computeFrictionReturnMagnitude(
        returnedMagnitude, unusedTrialDerivative, unusedBoundDerivative,
        state.tangentialNorm, state.frictionBound);
      state.tangentialTraction *= returnedMagnitude / state.tangentialNorm;
      if (frictionHardening > 0.0) {
        state.accumulatedPlasticSlip += std::max(
          0.0, (state.tangentialNorm - returnedMagnitude)
            / tangentialPenalty);
      }
    }
  }

  return state;
}


void
StructuralPenaltyContactBoundaryCondition::initializeFrom(const std::shared_ptr<InputRecord> &ir)
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
    // `friction`, `frictionTransition` and `frictionHardening` select the
    // experimental, development-phase friction branch (see the @warning on
    // the class declaration in the header) — leave `friction 0` for
    // production analyses until it has been fully verified.
    IR_GIVE_FIELD(ir, this->friction, _IFT_StructuralPenaltyContactBoundaryCondition_friction);
    IR_GIVE_OPTIONAL_FIELD(
      ir, this->frictionTransition,
      _IFT_StructuralPenaltyContactBoundaryCondition_frictionTransition);
    if (this->frictionTransition < 0.0 || this->frictionTransition >= 1.0) {
      OOFEM_ERROR(
        "StructuralPenaltyContactBoundaryCondition: frictiontransition must "
        "be in the interval [0, 1)");
    }
    IR_GIVE_OPTIONAL_FIELD(
      ir, this->frictionHardening,
      _IFT_StructuralPenaltyContactBoundaryCondition_frictionHardening);
    if (this->frictionHardening < 0.0 || this->frictionHardening >= 1.0) {
      OOFEM_ERROR(
        "StructuralPenaltyContactBoundaryCondition: frictionhardening must "
        "be in the interval [0, 1)");
    }
    if (this->frictionHardening > 0.0 && this->penalty_tangential <= 0.0) {
      OOFEM_ERROR(
        "StructuralPenaltyContactBoundaryCondition: frictionhardening "
        "requires a positive tangential penalty");
    }
    int nsd;
    IR_GIVE_OPTIONAL_FIELD(ir, nsd, _IFT_StructuralPenaltyContactBoundaryCondition_nsd);
    surface_dimension = nsd - 1;
    //
    int algo = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, algo, _IFT_StructuralPenaltyContactBoundaryCondition_algo);
    this->algo = algo;
    IR_GIVE_OPTIONAL_FIELD(ir, this->searchPadding, _IFT_StructuralPenaltyContactBoundaryCondition_searchPadding);
    if (this->searchPadding != -1.0 && this->searchPadding <= 0.0) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: searchpadding must be positive or omitted");
    }
    IR_GIVE_OPTIONAL_FIELD(ir, this->searchTol, _IFT_StructuralPenaltyContactBoundaryCondition_searchTol);
    if (this->searchTol < 0.0) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: searchtol must be nonnegative");
    }
    IR_GIVE_OPTIONAL_FIELD(ir, this->facetOwnershipHysteresis, _IFT_StructuralPenaltyContactBoundaryCondition_facetOwnershipHysteresis);
    if (this->facetOwnershipHysteresis < 0.0) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: facethysteresis must be nonnegative");
    }
    int featureFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(
      ir, featureFlag,
      _IFT_StructuralPenaltyContactBoundaryCondition_generalizedFeatures);
    this->generalizedFeatures = featureFlag != 0;
    int directionalProjectionFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(
      ir, directionalProjectionFlag,
      _IFT_StructuralPenaltyContactBoundaryCondition_directionalProjection);
    this->directionalProjection = directionalProjectionFlag != 0;
    if (this->directionalProjection && surface_dimension != 2) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: directionalprojection currently requires nsd 3");
    }
    if (this->directionalProjection && this->generalizedFeatures) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: directionalprojection and generalizedfeatures cannot be combined");
    }
    int autoPenaltyFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(
      ir, autoPenaltyFlag,
      _IFT_StructuralPenaltyContactBoundaryCondition_autoPenalty);
    this->autoPenalty = autoPenaltyFlag != 0;
    IR_GIVE_OPTIONAL_FIELD(
      ir, this->tangentMode,
      _IFT_StructuralPenaltyContactBoundaryCondition_tangentMode);
    if (this->tangentMode < 0 || this->tangentMode > 3) {
      OOFEM_ERROR(
        "StructuralPenaltyContactBoundaryCondition: tangentmode must be "
        "0 (automatic), 1 (rate-form analytical), 2 (finite difference), "
        "or 3 (finite-step analytical)");
    }
    if ((this->frictionTransition > 0.0 || this->frictionHardening > 0.0)
        && this->tangentMode == 1) {
      OOFEM_ERROR(
        "StructuralPenaltyContactBoundaryCondition: frictiontransition and "
        "frictionhardening require tangentmode 0, 2, or 3");
    }
    int fdCheck = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, fdCheck, _IFT_StructuralPenaltyContactBoundaryCondition_fdCheck);
    this->finiteDifferenceCheck = (fdCheck != 0);
    IR_GIVE_OPTIONAL_FIELD(ir, this->finiteDifferencePerturbation, _IFT_StructuralPenaltyContactBoundaryCondition_fdPerturbation);
    IR_GIVE_OPTIONAL_FIELD(ir, this->finiteDifferenceOutputPrefix, _IFT_StructuralPenaltyContactBoundaryCondition_fdOutputPrefix);
    IR_GIVE_OPTIONAL_FIELD(ir, this->finiteDifferenceRelativeTolerance, _IFT_StructuralPenaltyContactBoundaryCondition_fdTolerance);
}


void
StructuralPenaltyContactBoundaryCondition :: postInitialize()
{
  // take nodes of node set after everything is initialized
  this->masterSurfaceElements = domain->giveSet(this->masterSurfaceNumber)->giveElementList();
  this->slaveSurfaceElements = domain->giveSet(this->slaveSurfaceNumber)->giveElementList();
  ContactBoundaryCondition :: postInitialize();
  this->initializeAutomaticPenaltyFactors();
}

void
StructuralPenaltyContactBoundaryCondition :: saveContext(DataStream &stream, ContextMode mode)
{
  ContactBoundaryCondition :: saveContext(stream, mode);
  if (mode & CM_Definition) {
    if (!stream.write(penalty_normal)
        || !stream.write(penalty_tangential)
        || !stream.write(friction)
        || !stream.write(frictionTransition)
        || !stream.write(frictionHardening)
        || !stream.write(searchPadding)
        || !stream.write(searchTol)
        || !stream.write(facetOwnershipHysteresis)
        || !stream.write(generalizedFeatures)
        || !stream.write(directionalProjection)
        || !stream.write(autoPenalty)
        || !stream.write(tangentMode)
        || !stream.write(surface_dimension)
        || !stream.write(algo)
        || !stream.write(masterSurfaceNumber)
        || !stream.write(slaveSurfaceNumber)
        || !stream.write(finiteDifferenceCheck)
        || !stream.write(finiteDifferencePerturbation)
        || !stream.write(finiteDifferenceRelativeTolerance)
        || !stream.write(finiteDifferenceOutputPrefix)) {
      THROW_CIOERR(CIO_IOERR);
    }
  }
  if (mode & CM_State) {
    if (!contactSearchAlgorithm) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: contact search is not initialized while saving context");
    }
    contactSearchAlgorithm->saveContext(stream);
  }
}

void
StructuralPenaltyContactBoundaryCondition :: restoreContext(DataStream &stream, ContextMode mode)
{
  ContactBoundaryCondition :: restoreContext(stream, mode);
  if (mode & CM_Definition) {
    if (!stream.read(penalty_normal)
        || !stream.read(penalty_tangential)
        || !stream.read(friction)
        || !stream.read(frictionTransition)
        || !stream.read(frictionHardening)
        || !stream.read(searchPadding)
        || !stream.read(searchTol)
        || !stream.read(facetOwnershipHysteresis)
        || !stream.read(generalizedFeatures)
        || !stream.read(directionalProjection)
        || !stream.read(autoPenalty)
        || !stream.read(tangentMode)
        || !stream.read(surface_dimension)
        || !stream.read(algo)
        || !stream.read(masterSurfaceNumber)
        || !stream.read(slaveSurfaceNumber)
        || !stream.read(finiteDifferenceCheck)
        || !stream.read(finiteDifferencePerturbation)
        || !stream.read(finiteDifferenceRelativeTolerance)
        || !stream.read(finiteDifferenceOutputPrefix)) {
      THROW_CIOERR(CIO_IOERR);
    }

    slaveContactSurface = dynamic_cast<StructuralFEContactSurface *>(
      domain->giveContactSurface(slaveSurfaceNumber));
    masterContactSurface = dynamic_cast<StructuralFEContactSurface *>(
      domain->giveContactSurface(masterSurfaceNumber));
    if (!slaveContactSurface || !masterContactSurface) {
      OOFEM_ERROR("StructuralPenaltyContactBoundaryCondition: restored contact surfaces have incompatible types");
    }
    masterSurfaceElements = domain->giveSet(masterSurfaceNumber)->giveElementList();
    slaveSurfaceElements = domain->giveSet(slaveSurfaceNumber)->giveElementList();
  }

  if (!contactSearchAlgorithm) {
    this->setupContactSearchAlgorithm();
    contactSearchAlgorithm->createContactPairs();
  }
  if (autoPenalty && automaticPenaltyFactors.empty()) {
    this->initializeAutomaticPenaltyFactors();
  }
  if (mode & CM_State) {
    contactSearchAlgorithm->restoreContext(stream);
    // Context stores committed history, while the current master projection is
    // iteration-local. Rebuild it at the restored configuration so predictor
    // forces assembled before the first resumed nonlinear iteration include
    // the contact contribution from the checkpoint state.
    contactSearchAlgorithm->updateContactPairs(
      domain->giveEngngModel()->giveCurrentStep());
  }
}



void
StructuralPenaltyContactBoundaryCondition :: giveLocationArrays(std::vector< IntArray > &rows, std::vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
  rows.clear();
  cols.clear();

  // The sparse profile is built before the active master projection is known.
  // Reserve all admissible master/slave blocks so later contact activation or
  // facet switching cannot assemble outside a fixed Skyline profile.
  IntArray masterRow, masterColumn, slaveRow, slaveColumn;
  // Profile couplings depend on the slave facet, not on the number of
  // integration points on that facet.  Iterating contact pairs would emit the
  // same four blocks once per Gauss point (four times for a 2 x 2 rule), which
  // dominates startup and memory for large contact meshes.
  for (int slaveIndex = 1;
       slaveIndex <= slaveContactSurface->giveNumberOfContactElements();
       ++slaveIndex) {
    ContactElement *slaveElement =
      slaveContactSurface->giveContactElement_InSet(slaveIndex);
    slaveElement->giveLocationArray(slaveRow, dofs, r_s);
    slaveElement->giveLocationArray(slaveColumn, dofs, c_s);

    for (int masterIndex = 1; masterIndex <= masterContactSurface->giveNumberOfContactElements(); ++masterIndex) {
      ContactElement *masterElement = masterContactSurface->giveContactElement_InSet(masterIndex);
      masterElement->giveLocationArray(masterRow, dofs, r_s);
      masterElement->giveLocationArray(masterColumn, dofs, c_s);

      rows.push_back(masterRow);
      cols.push_back(masterColumn);
      rows.push_back(masterRow);
      cols.push_back(slaveColumn);
      rows.push_back(slaveRow);
      cols.push_back(masterColumn);
      rows.push_back(slaveRow);
      cols.push_back(slaveColumn);
    }
  }

  // A convected friction state crossing a facet edge depends on both the new
  // and committed master facets. Reserve the rectangular master-master blocks
  // for topologically adjacent facets. Slave-old-master blocks are already
  // covered by the all-slave/all-master loop above.
  const int nMasterElements = masterContactSurface->giveNumberOfContactElements();
  for (int i = 1; i <= nMasterElements; ++i) {
    ContactElement *first = masterContactSurface->giveContactElement_InSet(i);
    for (int j = i + 1; j <= nMasterElements; ++j) {
      ContactElement *second = masterContactSurface->giveContactElement_InSet(j);
      int sharedNodes = 0;
      for (int firstNode = 1; firstNode <= first->giveNumberOfNodes(); ++firstNode) {
        for (int secondNode = 1; secondNode <= second->giveNumberOfNodes(); ++secondNode) {
          if (first->giveNode(firstNode) == second->giveNode(secondNode)) {
            ++sharedNodes;
            break;
          }
        }
      }
      if (sharedNodes < surface_dimension) {
        continue;
      }

      first->giveLocationArray(masterRow, dofs, r_s);
      second->giveLocationArray(masterColumn, dofs, c_s);
      rows.push_back(masterRow);
      cols.push_back(masterColumn);

      second->giveLocationArray(masterRow, dofs, r_s);
      first->giveLocationArray(masterColumn, dofs, c_s);
      rows.push_back(masterRow);
      cols.push_back(masterColumn);
    }
  }
}


  

} // namespace oofem
