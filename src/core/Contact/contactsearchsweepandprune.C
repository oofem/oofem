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

#include "contactsearchsweepandprune.h"
#include "floatarrayf.h"
#include <algorithm>

namespace oofem {

namespace {

constexpr bool
boundLess(const Bound &a, const Bound &b)
{
	if (a.value < b.value) {
		return true;
	}
	if (b.value < a.value) {
		return false;
	}
	if (a.isMin != b.isMin) {
		return a.isMin; // Closed intervals overlap when their endpoints touch.
	}
	return a.id < b.id;
}

static_assert(boundLess(Bound{0.0, 2u, false, true}, Bound{0.0, 1u, true, false}),
              "SAP touching intervals require minimum endpoints before maximum endpoints");
static_assert(!boundLess(Bound{0.0, 1u, true, false}, Bound{0.0, 2u, false, true}),
              "SAP touching endpoint ordering must be antisymmetric");
static_assert(boundLess(Bound{0.0, 1u, false, true}, Bound{0.0, 2u, true, true}),
              "SAP equal endpoint ties require deterministic entity ordering");

}

ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune :: ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune(
	FEContactSurface *scs,
	FEContactSurface *mcs,
	Domain *d
) : ContactSearchAlgorithm_Surface2FESurface_3d(scs, mcs, d)
{
	isInitialized = false;
}

void
ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune :: initialize(TimeStep *tStep)
{
	unsigned int n_slaves = contactPairs.size();
	unsigned int n_masters = this->masterContactSurface->giveNumberOfContactElements();
	unsigned int n = n_slaves + n_masters;
	//
	aabbs.resize(n_slaves + n_masters);
	updateAABBs(tStep);
	//
	this->potentialPairs.clear();
	for (int axis=0; axis<3; axis++) {
		this->boundss[axis] = {};
		this->potentialPairsPerAxes[axis].clear();
	}
	for (unsigned int i = 0; i < n; ++i) {
		bool isSlave = i >= n_masters;
		const auto& aabb = aabbs[i];
		for (unsigned int axis = 0; axis < 3; ++axis) {
			auto& bounds = this->boundss[axis];
			bounds.push_back({
				aabb.min[axis],
				i,
				isSlave,
				true
			});
			bounds.push_back({
				aabb.max[axis],
				i,
				isSlave,
				false
			});
		}
	}
	for (int axis = 0; axis < 3; ++axis) {
		auto& bounds = this->boundss[axis];
		auto& potentialPairsPerAxis = this->potentialPairsPerAxes[axis];
		potentialPairsPerAxis.clear();
		std::sort(bounds.begin(), bounds.end(), boundLess);
		const size_t nBounds = bounds.size();
		for (size_t i = 0; i + 1 < nBounds; ++i) {
			const auto& bi = bounds[i];
			if (!bi.isMin) continue;
			for (size_t j = i + 1; j < nBounds; ++j) {
				const auto& bj = bounds[j];
				//
				if (bi.id == bj.id) break;
				//
				if (!bj.isMin) continue;
				if (bi.isSlave == bj.isSlave) continue;
				//
				unsigned int i_master;
				unsigned int i_slave;
				if (bi.isSlave) {
					i_slave = bi.id;
					i_master = bj.id;
				} else {
					i_slave = bj.id;
					i_master = bi.id;
				}
				IJ ij = {i_master, i_slave};
				//
				potentialPairsPerAxis.insert(ij);
			}
		}
	}
	//
	const auto& potentialPairsPerAxis0 = this->potentialPairsPerAxes[0];
	const auto& potentialPairsPerAxis1 = this->potentialPairsPerAxes[1];
	const auto& potentialPairsPerAxis2 = this->potentialPairsPerAxes[2];
	for (const auto& ij : potentialPairsPerAxis0) {
		if (
			contains(potentialPairsPerAxis1, ij)
			&&
			contains(potentialPairsPerAxis2, ij)
		) {
			this->addPotentialPair(ij);
		}
	}
	this->isInitialized = true;
}

void
ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune :: updateAABBs(TimeStep *tStep)
{
	unsigned int index = 0;
	double searchPadding = this->computeSearchPadding(tStep);
	//
	unsigned int n_masters = this->masterContactSurface->giveNumberOfContactElements();
	for (unsigned int i = 1; i <= n_masters; ++i) {
		auto aabb = this->computeSweptMasterSearchAABB(i, tStep);
		aabb.expandBy(searchPadding);
		aabbs[index] = aabb;
		index++;
	}
	//
	for (const auto& slave : contactPairs) {
		auto aabb = slave->computeSweptSlaveAABB(tStep);
		aabb.expandBy(searchPadding);
		aabbs[index] = aabb;
		index++;
	}

}

void
ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune :: updateBoundss(TimeStep *tStep)
{
	updateAABBs(tStep);
	//
	for (int axis = 0; axis < 3; ++axis) {
		auto& bounds =
			(axis == 0) ? std::get<0>(this->boundss)
			:
			(axis == 1) ? std::get<1>(this->boundss)
			:
			std::get<2>(this->boundss)
			;
		for (auto& b : bounds) {
			const auto& aabb = aabbs[b.id];
			b.value = b.isMin ? aabb.min[axis] : aabb.max[axis];
		}
	}
}

void
ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune :: insertionSort()
{
	for (int axis = 0; axis < 3; ++axis) {
		auto& bounds = this->boundss[axis];
		for (size_t i = 1; i < bounds.size(); ++i) {
			for (size_t j = i; j >= 1; --j) {
				Bound& b1 = bounds[j - 1];
				Bound& b2 = bounds[j];
				if (!boundLess(b2, b1)) break;
				std::swap(b1, b2);
				this->solveInversion(axis, b1, b2);
			}
		}
	}
}

void
ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune :: addPotentialPair(IJ ij)
{
	this->potentialPairs.insert(ij);
}

void
ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune :: deletePotentialPair(IJ ij)
{
	this->potentialPairs.erase(ij);
}

void
ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune :: solveInversion(
	int axis,
	const Bound& b1,
	const Bound& b2
)
{
	if (b1.isSlave == b2.isSlave) return;
	if (b1.isMin == b2.isMin) return;
	int axis2, axis3;
	if (axis == 0) {
		axis2 = 1;
		axis3 = 2;
	} else if (axis == 1) {
		axis2 = 0;
		axis3 = 2;
	} else {
		axis2 = 0;
		axis3 = 1;
	}
	auto& potentialPairsPerAxis1 = this->potentialPairsPerAxes[axis];
	auto& potentialPairsPerAxis2 = this->potentialPairsPerAxes[axis2];
	auto& potentialPairsPerAxis3 = this->potentialPairsPerAxes[axis3];
	unsigned int i_master;
	unsigned int i_slave;
	if (b1.isSlave) {
		i_slave = b1.id;
		i_master = b2.id;
	} else {
		i_slave = b2.id;
		i_master = b1.id;
	}
	IJ ij = {i_master, i_slave};
	if (b1.isMin) {
		potentialPairsPerAxis1.insert(ij);
		if (
			contains(potentialPairsPerAxis2, ij)
			&&
			contains(potentialPairsPerAxis3, ij)
		) {
			this->addPotentialPair(ij);
		}
	} else {
		potentialPairsPerAxis1.erase(ij);
		this->deletePotentialPair(ij);
	}
}

void
ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune :: updateContactPairs(TimeStep *tStep)
{
	if (!this->isInitialized) {
		this->initialize(tStep);
	} else {
		this->updateBoundss(tStep);
		this->insertionSort();
	}
	//
	std::map<unsigned int, std::vector<unsigned int>> slave2masters;
	for (IJ ij : potentialPairs) {
		auto [i_master, i_slave] = ij;
		if (slave2masters.find(i_slave) == slave2masters.end()) {
			slave2masters[i_slave] = {};
		}
		slave2masters[i_slave].push_back(i_master);
	}
	//
	const unsigned int n_masters = this->masterContactSurface->giveNumberOfContactElements();
	//
		for (size_t i_slave = 0; i_slave < contactPairs.size(); ++i_slave) {
			auto& cp = contactPairs[i_slave];
			std::vector<int> masterCandidateIndices;

		const unsigned int encodedSlaveId = n_masters + static_cast<unsigned int>(i_slave);
		auto candidateIt = slave2masters.find(encodedSlaveId);
		if (candidateIt != slave2masters.end()) {
			for (unsigned int i_master_0 : candidateIt->second) {
				if (i_master_0 < n_masters) {
					masterCandidateIndices.push_back(i_master_0 + 1);
				}
			}
		}

		this->appendCurrentActiveMasterCandidate(cp.get(), masterCandidateIndices);
			if (masterCandidateIndices.empty()) {
				cp->clearContactState();
				continue;
			}

			this->updatePairFrom3dCandidates(cp.get(), masterCandidateIndices, tStep);
		}
}

void
ContactSearchAlgorithm_Surface2FESurface_3d_SweepAndPrune :: restoreContext(DataStream &stream)
{
	ContactSearchAlgorithm_Surface2FESurface_3d :: restoreContext(stream);
	potentialPairs.clear();
	for (auto &pairs : potentialPairsPerAxes) {
		pairs.clear();
	}
	for (auto &bounds : boundss) {
		bounds.clear();
	}
	aabbs.clear();
	isInitialized = false;
}

} // end namespace oofem
