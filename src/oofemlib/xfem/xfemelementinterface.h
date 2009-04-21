#ifndef xfemelementinterface_h
#define xfemelementinterface_h

#include "interface.h"
#include "gausspnt.h"
#include "xfemmanager.h"
#include "matresponsemode.h"
#include "structuralelement.h"

/** provides xfem interface for an element */
class XfemElementInterface : public Interface
{
public:
    /// constructor
    XfemElementInterface(Element *e) : Interface() { 
        this->element = e; }
    /// creates enriched part of B matrix
    void XfemElementInterface_createEnrBmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    /// partitions the element into patches by a triangulation
    void XfemElementInterface_partitionElement(AList< Triangle > *answer, AList< FloatArray > *together);
    /// updates integration rule based on the triangulation
    //void XfemElementInterface_updateIntegrationRule();
    /// helpful routine to put the nodes for triangulation together, should be in protected members probably
    void XfemElementInterface_prepareNodesForDelaunay(AList< FloatArray > *answer1, AList< FloatArray > *answer2);
    /*virtual void XfemElementInterface_computeBmatrixAt(GaussPoint *, FloatMatrix &answer,
                                   int lowerIndx = 1, int upperIndx = ALL_STRAINS) = 0;
    virtual int XfemElementInterface_computeNumberOfDofs(EquationID ut) = 0;
    virtual void XfemElementInterface_giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const = 0;
    virtual void XfemElementInterface_computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *, TimeStep *tStep) = 0;
    virtual double XfemElementInterface_computeVolumeAround(GaussPoint *) = 0;
    virtual void XfemElementInterface_computeVectorOf(EquationID type, ValueModeType u, TimeStep *stepN, FloatArray &answer) = 0;
    virtual void XfemElementInterface_giveLocationArray(IntArray & locationArray, EquationID) const = 0;
    void computeBmatrixAt(GaussPoint *, FloatMatrix &answer,
                                   int lowerIndx = 1, int upperIndx = ALL_STRAINS);
    int computeNumberOfDofs(EquationID ut);
    void giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const ;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *, TimeStep *tStep) ;
    double computeVolumeAround(GaussPoint *);
    void computeVectorOf(EquationID type, ValueModeType u, TimeStep *stepN, FloatArray &answer) ;
    void giveLocationArray(IntArray & locationArray, EquationID) const; */
    void XfemElementInterface_updateIntegrationRule();
protected:
    Element *element;
};
#endif
