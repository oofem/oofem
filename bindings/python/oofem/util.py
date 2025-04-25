try: from . import oofempy         # when installed
except ImportError: import oofempy # imported in the build-tree, oofempy being in PYTHONPATH
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

def plot2Dmesh(d, tstep, field=None, fieldValueIndex=0):
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
            z.append(val[fieldValueIndex])

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


def plot1Dmesh(d, tstep, xind=0, yind=1, evals=None, warpField=None, warpScale=1.0, label="", nodeLabels=False, elementLabels=False):
    x=[];
    y=[];
    
    fig, ax = plt.subplots()
    #plt.set_title('3D OOFEM mesh')
    ax.set_aspect('equal')

    if (evals):
        # Normalize values for color mapping
        mmin = min(evals)
        mmax = max(evals)
        norm = plt.Normalize(mmin, mmax)
        colors = plt.cm.viridis(norm(evals))        

    ElemBbox = dict(boxstyle='round', fc='lavender', ec='blue', alpha=0.5)
    NodeBbox = dict(boxstyle='round', fc='mistyrose', ec='red', alpha=0.5)

    for i in range(1, d.giveNumberOfElements()+1):
        e=d.giveElement(i)
        gt=e.giveGeometryType()
        if (gt == oofempy.Element_Geometry_Type.EGT_line_1):
            d1 = e.giveDofManager(1)
            d2 = e.giveDofManager(2)
            c1 = d1.giveCoordinates()
            c2 = d2.giveCoordinates()
            
            
            if (warpField):
                val1 = oofempy.FloatArray()
                warpField.evaluateAt(val1, d1, oofempy.ValueModeType.VM_Total, tstep) 
                val2 = oofempy.FloatArray()
                warpField.evaluateAt(val2, d2, oofempy.ValueModeType.VM_Total, tstep) 
                x = (c1[xind] + val1[xind]*warpScale, c2[xind] + val2[xind]*warpScale)
                y = (c1[yind] + val1[yind]*warpScale, c2[yind] + val2[yind]*warpScale)

            else:
                x  = (c1[xind], c2[xind])
                y  = (c1[yind], c2[yind])

            if (evals):
                plt.plot(x,y, color=colors[i-1], marker='s', linestyle='solid', linewidth=2, markersize=4, mfc='red')
            else:
                plt.plot(x,y, color='black', marker='s', linestyle='solid', linewidth=2, markersize=4, mfc='red')
            if (elementLabels):
                plt.text((x[0]+x[1])/2, (y[0]+y[1])/2, str(e.giveNumber()), fontsize=8, ha='center', va='center', bbox=ElemBbox)

    if (nodeLabels):
        for i in range(1, d.giveNumberOfDofManagers()+1):
            dm = d.giveDofManager(i)
            coords = dm.giveCoordinates()
            if (warpField):
                val = oofempy.FloatArray()
                warpField.evaluateAt(val, dm, oofempy.ValueModeType.VM_Total, tstep) 
                coords[xind] += val[xind]*warpScale
                coords[yind] += val[yind]*warpScale
            plt.text(coords[xind], coords[yind], str(dm.giveNumber()), fontsize=8, ha='center', va='center', bbox=NodeBbox)

    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, label=label, ax=ax, orientation='horizontal')   
    return plt
