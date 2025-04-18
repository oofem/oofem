import oofempy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


#
# Initializes domain with given sequences of nodes, elements, .....
#
def setupDomain (domain, nodes, elems, css, mats, bcs, ics, ltfs, sets):
    
    domain.resizeDofManagers(len(nodes)) # process nodes
    for n in nodes:
        domain.setDofManager(n.number, n)
    
    domain.resizeElements(len(elems)) # process elements
    for e in elems:
        domain.setElement(e.number, e)
    
    domain.resizeCrossSectionModels(len(css)) #process cross sections
    for c in css:
        domain.setCrossSection(c.number, c)
    
    domain.resizeMaterials(len(mats)) # process materials
    for m in mats:
        domain.setMaterial(m.number, m)
    
    domain.resizeBoundaryConditions(len(bcs)) #process boundary conditions
    for bc in bcs:
        domain.setBoundaryCondition(bc.number, bc)
    
    domain.resizeInitialConditions(len(ics)) #process initial conditions
    for ic in ics:
        domain.setInitialCondition(ic.number, ic)
    
    domain.resizeFunctions(len(ltfs)) # process (load time) functions
    for ltf in ltfs:
        domain.setFunction(ltf.number, ltf)
    
    domain.resizeSets(len(sets)) # process sets
    for s in sets:
        domain.setSet(s.number, s)
       
    return

def plot2D(d, field=None, fieldValueIndex=0, tstep=None):
    x=[];
    y=[];
    z=[];
    elements=[];
    for i in range(1,d.giveNumberOfDofManagers()+1):
        coords = d.giveDofManager(i).giveCoordinates()
        coords.resize(2)
        x.append(coords[0]);
        y.append(coords[1]);
        if (field):
            val = oofempy.FloatArray(6)
            field.evaluateAt(val, coords, oofempy.ValueModeType.VM_Total, tstep) 
            z.append(val[0])

    #print(nodes)
    for i in range(1, d.giveNumberOfElements()+1):
        e=d.giveElement(i)
        gt=e.giveGeometryType()
        if (gt == oofempy.Element_Geometry_Type.EGT_triangle_1):
            elements.append((e.giveDofManager(1).giveNumber()-1, e.giveDofManager(2).giveNumber()-1, e.giveDofManager(3).giveNumber()-1))
    triangulation = tri.Triangulation(x, y, elements)

    fig, ax = plt.subplots()
    #plt.set_title('3D OOFEM mesh')
    ax.set_aspect('equal')

    if (field):
        levels = np.linspace(min(z), max(z), 12)
        plt.tricontourf(triangulation, z, levels=levels, alpha=1.0)
        plt.colorbar(orientation='horizontal')
    else:
        plt.triplot(triangulation, linewidth=0.2, color='black', alpha=1.0)
    
    return plt