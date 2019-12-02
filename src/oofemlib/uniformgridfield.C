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

#include <cassert>
#include <cmath>
#include "uniformgridfield.h"
#include "dofmanager.h"
#include "error.h"

#include <iostream>

namespace oofem {

// inspired by https://github.com/woodem/woo/blob/master/pkg/dem/FlowAnalysis.hpp#L20
// returns lo-indices of containing cell (clamped to domain)
// and also computes normalized (0…1)×..×(0…1) coords within that cell 
void UniformGridField::xyz2ijk(const FloatArray& xyz, IntArray& ijk, FloatArray& normXyz) const {
    assert(xyz.giveSize()<=lo.giveSize());
    assert(hi.giveSize()==lo.giveSize());
    assert(div.giveSize()==lo.giveSize());
    ijk.resize(xyz.giveSize());
    normXyz.resize(xyz.giveSize());
    for(int ax=0; ax<xyz.giveSize(); ax++){
        // FIXME: cellDim should be cached and not computed over and over
        double cellDim=(hi[ax]-lo[ax])/div[ax];
        double t=(xyz[ax]-lo[ax])/cellDim;
        // clamp to grid size
        if(t>=div[ax]){ ijk[ax]=div[ax]-1; normXyz[ax]=1.; }
        else if(t<0) { ijk[ax]=0;          normXyz[ax]=0.; }
        else         { ijk[ax]=(int)std::floor(t); normXyz[ax]=t-ijk[ax]; }
    }
}

void UniformGridField::setGeometry(const FloatArray& lo_, const FloatArray& hi_, const IntArray& div_){
    // std::cerr<<"lo dimension: "<<lo.giveSize()<<std::endl;
    if(lo_.giveSize()!=hi_.giveSize()) OOFEM_ERROR("lo and hi must have identical dimension.");
    if(lo_.giveSize()!=div_.giveSize()) OOFEM_ERROR("lo and div must have identical dimension.");
    if(lo_.giveSize()!=2 && lo_.giveSize()!=3) OOFEM_ERROR("Dimension must be 2 or 3.");
    for(int i=0; i<lo_.giveSize(); i++){
        if(lo_[i]>=hi_[i]) OOFEM_ERROR("lo must be strictly smaller than hi along all axes.");
        if(div_[i]<1) OOFEM_ERROR("div must be at least 1 along all axes.")
    }
    lo=lo_;
    hi=hi_;
    div=div_;
    #if 0
        this->precomputeInternal();
    #endif
}
void UniformGridField::setValues(const FloatArray& vv){
    int s=1;
    for(int i=0; i<div.giveSize(); i++) s*=div[i]+1;
    if(vv.giveSize()!=s) OOFEM_ERROR((std::string("Array size must be exactly prod(div[i]+1)=")+std::to_string(s)).c_str());
    values=vv;
}

double UniformGridField::nodeValue2d(int i, int j){
    assert(div.giveSize()==2);
    assert(values.giveSize()==(div[0]+1)*(div[1]+1));
    return values[(div[1]+1)*i+j];
}

double UniformGridField::nodeValue3d(int i, int j, int k){
    assert(div.giveSize()==3);
    assert(values.giveSize()==(div[0]+1)*(div[1]+1)*(div[2]+1));
    return values[(div[2]+1)*(div[1]+1)*i+(div[2]+1)*j+k];
}

// see https://github.com/woodem/woo/blob/master/pkg/dem/FlowAnalysis.cpp
// and https://woodem.org/user/flow-analysis.html for explanation about the interpolation routine
int UniformGridField::evaluateAt(FloatArray &answer, const FloatArray &coords,
                           ValueModeType mode, TimeStep *tStep){
    // scalar value
    answer.resize(1);
    double& ret(answer[0]);
    ret=0;

    // find cell containing coords, and coords within the cell
    IntArray ijk; FloatArray normXyz;
    this->xyz2ijk(coords,ijk,normXyz);

    // 3D interpolation
    if(ijk.giveSize()==3){
        const int& i(ijk[0]); const int& j(ijk[1]); const int& k(ijk[2]);
        int I(i+1), J(j+1), K(k+1);
        const double& x(normXyz[0]); const double& y(normXyz[1]); const double& z(normXyz[2]);
        double X(1-x), Y(1-y), Z(1-z);
        double weights[]={
            X*Y*Z,x*Y*Z,x*y*Z,X*y*Z,
            X*Y*z,x*Y*z,x*y*z,X*y*z
        };
        int pts[][3]={
            {i,j,k},{I,j,k},{I,J,k},{i,J,k},
            {i,j,K},{I,j,K},{I,J,K},{i,J,K}
        };
        assert(abs(weights[0]+weights[1]+weights[2]+weights[3]+weights[4]+weights[5]+weights[6]+weights[7]-1) < 1e-5);
        for(int p=0; p<8; p++){
            // std::cerr<<p<<": at "<<coords[0]<<","<<coords[1]<<","<<coords[2]<<", ijk "<<i<<","<<j<<","<<k<<" normXYZ "<<x<<","<<y<<","<<z<<std::endl;
            ret+=weights[p]*this->nodeValue3d(pts[p][0],pts[p][1],pts[p][2]);
        }
    }
    // 2D interpolation
    else if (ijk.giveSize()==2){
        const int& i(ijk[0]); const int& j(ijk[1]); int I(i+1); int J(j+1);
        const double& x(normXyz[0]); const double y(normXyz[1]); double X(1-x); double Y(1-y);
        double weights[]={X*Y,x*Y,X*y,x*y}; // TODO: check
        int pts[][2]={{i,j},{I,j},{i,J},{I,J}}; // TODO: check
        assert(abs(weights[0]+weights[1]+weights[2]+weights[3]-1)<1e-5);
        for(int p=0; p<4; p++){
            // std::cerr<<p<<": at "<<coords[0]<<","<<coords[1]<<" ij "<<i<<","<<j<<" normXY "<<x<<","<<y<<std::endl;
            ret+=weights[p]*this->nodeValue2d(pts[p][0],pts[p][1]); 
        }
    }
    // any other is unsupported
    else {
        OOFEM_ERROR((std::string("UniformGridField::evaluateAt: erroneous dimension of input coordinates (")+std::to_string(coords.giveSize())+", must be 2 or 3).").c_str());
        return 1;
    }

    return 0; // OK
}



int
UniformGridField :: evaluateAt(FloatArray &answer, DofManager *dman, ValueModeType mode, TimeStep *tStep)
{
    return evaluateAt(answer, dman->giveCoordinates(), mode, tStep) == 1;
}



} // end namespace oofem
