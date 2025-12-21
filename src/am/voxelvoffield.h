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

#ifndef fvoxelvoffield_h
#define fvoxelvoffield_h

#include "field.h"
#include "element.h"
#include "node.h"
#include "VoxelGrid.h"
#include "timestep.h"
#include "floatarray.h"


namespace oofem {
class VoxelVOFField : public oofem::Field
{
    protected:
    VoxelGrid *voxelGrid;
    public:
    VoxelVOFField(VoxelGrid *vg = nullptr) : Field (oofem::FieldType::FT_VOF), voxelGrid(vg) {}
    void setGrid(VoxelGrid *vg) { voxelGrid = vg; }
    int evaluateAt(FloatArray &answer, const FloatArray &coords, ValueModeType mode, TimeStep *tStep) override {
        answer.resize(1);
        auto indices = voxelGrid->get_indices_from_point({coords[0], coords[1], coords[2]});
        int indx = voxelGrid->get_index( std::get<0>(indices),
                                     std::get<1>(indices),
                                     std::get<2>(indices) );
        if (voxelGrid->is_active(indx)) {
            answer[0]=voxelGrid->get_voxel(indx).getVofAt(tStep->giveIntrinsicTime());
        } else {
            answer[0] = 0.0; // empty voxel
        }
        return 0;
    }

    int evaluateAt(FloatArray &answer, DofManager *dman,
                   ValueModeType mode, TimeStep *tStep) override {
        answer.resize(1);
        answer[0]=0.0;
        return 0;
    }

    int evaluateAt(FloatArray &answer, Element *el,
                   ValueModeType mode, TimeStep *tStep) override {
        // there is no mapping from element number to voxel idx, so compute element center and get voxel from that
        FloatArray coords(3);
        int nnodes = el->giveNumberOfNodes();
        for (int i=1; i<=nnodes; i++) {
            auto c = el->giveNode(i)->giveCoordinates();
            coords.add(c);
        }
        coords *= 1000./nnodes; // convert from m to mm
        auto indices = voxelGrid->get_indices_from_point({coords[0], coords[1], coords[2]});
        int indx = voxelGrid->get_index( std::get<0>(indices),
                                     std::get<1>(indices),
                                     std::get<2>(indices) );
        if (voxelGrid->is_active(indx)) {
            answer[0]=voxelGrid->get_voxel(indx).getVofAt(tStep->giveIntrinsicTime());
        } else {
            answer[0] = 0.0; // empty voxel
        }
        return 0;
    }
   
    void saveContext(DataStream &stream) override {};
    void restoreContext(DataStream &stream) override {};


     /// @return Class name of the receiver.
    virtual const char *giveClassName() const override {
        return "VoxelVOFField";
    }

    // for Field classes supporting instantiation from input record
    virtual void initializeFrom(InputRecord &ir) override { };
};
} // end namespace oofem
#endif // fvoxelvoffield_h