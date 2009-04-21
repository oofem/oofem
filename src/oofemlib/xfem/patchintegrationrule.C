#include "patchintegrationrule.h"
#include "xfemelementinterface.h"
#include "patch.h"
#include "integrationrule.h"
#include "gaussintegrationrule.h"

PatchIntegrationRule :: PatchIntegrationRule(int n, Element *e, AList< Triangle > *triangles) : IntegrationRule(n, e) {
    this->patchSet = new AList< Patch >(0);
    for ( int i = 1; i  <= triangles->giveSize(); i++ ) {
        // too specific here, factory method outside
        Patch *patch = new TrianglePatch(elem);
        for ( int j = 1; j <= triangles->at(i)->giveVertices()->giveSize(); j++ ) {
            FloatArray *nCopy = new FloatArray( *triangles->at ( i )->giveVertex ( j ) );
            patch->setVertex(nCopy);
        }

        this->patchSet->put(i, patch);
    }
}

PatchIntegrationRule :: ~PatchIntegrationRule() {
    delete patchSet;
}

void PatchIntegrationRule :: computeGpsForPatch(Patch *patch, AList< GaussPoint > *gps, MaterialMode matMode) {
    GaussIntegrationRule tmpIR(1, elem, 1, 1);
    // too specific here, factory method outside
    tmpIR.setUpIntegrationPoints(_Triangle, patch->giveNrVertices(), matMode);
    for ( int j = 0; j < tmpIR.getNumberOfIntegrationPoints(); j++ ) {
        GaussPoint *gp = tmpIR.getIntegrationPoint(j);
        double weight = 0;
        FloatArray *coor = new FloatArray();
        GaussPoint *cp = new GaussPoint(this, 1, coor, weight, matMode);
        FloatArray *gpCoor = gp->giveCoordinates();
        cp->setCoordinates(* gpCoor);
        cp->setWeight( gp->giveWeight() );
        gps->put(j + 1, cp);
    }
}

void PatchIntegrationRule :: computeGps(MaterialMode matMode) {
    AList< GaussPoint >nGparray;
    int count = 0;
    for ( int i = 1; i <= patchSet->giveSize(); i++ ) {
        Patch *p = this->patchSet->at(i);
        AList< GaussPoint >gps;
        computeGpsForPatch(p, & gps, matMode);
        for ( int j = 1; j <= gps.giveSize(); j++ ) {
            count++;
            p->convertGPIntoParental( gps.at(j) );
            nGparray.put( count, gps.at(j) );
            gps.unlink(j);
        }
    }

    this->numberOfIntegrationPoints = count;
    this->gaussPointArray = new GaussPoint * [ numberOfIntegrationPoints ];
    for ( int i = 0; i < count; i++ ) {
        this->gaussPointArray [ i ] = nGparray.at(i + 1);
        nGparray.unlink(i + 1);
    }
}

int PatchIntegrationRule :: SetUpPointsOnSquare(int nPoints, MaterialMode matMode, GaussPoint ***gps) {
    this->computeGps(matMode);
    return nPoints;
}
