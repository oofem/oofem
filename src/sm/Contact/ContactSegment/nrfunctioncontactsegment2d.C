#include "nrfunctioncontactsegment2d.h"

namespace oofem {

    void NRFunctionContactSegment2D::computeContactPoint(FloatArray & answer, FloatArray& normal, const FloatArray & nodeCoords)
    {
        if ( nodeCoords.giveSize() != 2 ) {
            OOFEM_ERROR("An incompatible coordinate size (%i) encountered. Algorithm is only for 2D functions.", nodeCoords.giveSize());
        }

        double x_node = nodeCoords.at(1);
        double y_node = nodeCoords.at(2);

        //minimizing the distance function d(x) = sqrt((x - x_node)^2 + (f(x) - y_node)^2)
        //omit the unnecessary root and the x_node^2 and y_node^2 to get
        //d(x) = x^2 - 2*x*x_node + f(x)^2 - 2*f(x)*y_node
        //dd(x)/dx = 2*x - 2*x_node + 2* f(x) * df(x)/dx - 2*y_node*df(x)/dx
        //d^2d(x)/dx^2 = 2 + 2*f(x)*(d^2f(x)/dx^2) + 2*(df(x)/dx)^2 - 2*y_node*(d^2f(x)/dx^2)

        //Newton-Raphson, inital guess x_node
        double x_c = x_node;
        double g, h;
        int k = 0; //iterator
        int maxiter = NRFunctionContact_Maxiter;
        double tol = NRFunctionContact_Tolerance;
        while ( k < maxiter ) {
            g = 2. * x_c - 2. * x_node + 2. * functionValue(x_c)*derivativeValue(x_c) - 2. * y_node*derivativeValue(x_c);
            if ( abs(g) <= tol ) break;
            h = 2. + 2.*functionValue(x_c)*doubleDerivativeValue(x_c) + 2.*derivativeValue(x_c)*derivativeValue(x_c) - 2.*y_node*doubleDerivativeValue(x_c);
            x_c += -g / h;
            k++;
        }

        if ( k >= maxiter ) OOFEM_WARNING("Searching for contact with analytical function: Newton-Rhapson method did not converge in %i iterations. Continuing.", k);

        //we found the contact point
        //FloatArray contactPoint(2); -- now contact point is the answer and normal is an additional answer
        answer.resize(2);
        answer.at(1) = x_c;
        answer.at(2) = functionValue(x_c);

        //now we have contact point, what we need is the normal at the point of contact
        //can be obtained as a derivative of the function
        //the tangent is in the form
        // y = slope*x
        double slope = derivativeValue(x_c);
        //therefore
        // slope*x - y = 0
        //and (slope, -1) is the normal vector to the tangent line
        normal.resize(2);
        normal.at(1) = slope;
        normal.at(2) = -1.;
        normal.normalize();

        normal.times(-1.);
    }
}