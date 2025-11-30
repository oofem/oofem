try: from . import oofempy         # when installed
except ImportError: import oofempy # imported in the build-tree, oofempy being in PYTHONPATH

import numpy as np
import matplotlib
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
    
    domain.initializeFinish()
       
    return

def plot2Dmesh(ax, d, tstep, field=None, fieldValueIndex=0, warpField=None, warpScale=1.0):
    """ Produces a 2D mesh plot on the given Axes object 
    
    Parameters:
    ax : matplotlib.axes.Axes
        The Axes object on which to plot the 2D mesh.
    d : oofem.Domain
        The domain containing the mesh data.
    tstep : oofem.TimeStep
        The time step for which to plot the mesh.
    field : oofem.Field, optional
        The field to evaluate and plot on the mesh.
    fieldValueIndex : int, optional
        The index of the field value to plot.
    warpField : oofem.Field, optional
        The field used to warp the mesh coordinates.
    warpScale : float, optional
        The scale factor for warping the mesh coordinates.
    """
    x=[];
    y=[];
    z=[];
    elements=[];
    for i in range(1,d.giveNumberOfDofManagers()+1):
        dman = d.giveDofManager(i)
        coords = dman.giveCoordinates()
        coords.resize(2)
        
        if (field):
            val = oofempy.FloatArray(6)
            try:
                field.evaluateAt(val, dman, oofempy.ValueModeType.VM_Total, tstep) 
                z.append(val[fieldValueIndex])
            except:
                print("field querry failed")
                pass


        if (warpField):
            val = oofempy.FloatArray(6)
            warpField.evaluateAt(val, dman, oofempy.ValueModeType.VM_Total, tstep) 
            coords[0] += val[0]*warpScale
            coords[1] += val[1]*warpScale
        x.append(coords[0]);
        y.append(coords[1]);

        
    #print(nodes)
    for i in range(1, d.giveNumberOfElements()+1):
        e=d.giveElement(i)
        gt=e.giveGeometryType()
        if (gt == oofempy.Element_Geometry_Type.EGT_triangle_1):
            elements.append((e.giveDofManager(1).giveNumber()-1, e.giveDofManager(2).giveNumber()-1, e.giveDofManager(3).giveNumber()-1))
        elif (gt == oofempy.Element_Geometry_Type.EGT_quad_1):
            elements.append((e.giveDofManager(1).giveNumber()-1, e.giveDofManager(2).giveNumber()-1, e.giveDofManager(3).giveNumber()-1))
            elements.append((e.giveDofManager(1).giveNumber()-1, e.giveDofManager(3).giveNumber()-1, e.giveDofManager(4).giveNumber()-1))
            
    triangulation = tri.Triangulation(x, y, elements)
    
    if (field):
        levels = np.linspace(int(min(z)-1.0), int(max(z)+1.0), 20)
        plt.tricontourf(triangulation, z, levels=levels, alpha=1.0)
        plt.colorbar(orientation='horizontal')
    
    plt.triplot(triangulation, linewidth=0.2, color='black', alpha=0.5)
    
    return plt


def plot1Dmesh(ax, d, tstep, xind=0, yind=1, evals=None, warpField=None, warpScale=1.0, label="", nodeLabels=False, elementLabels=False):
    """ Produces a 1D mesh plot on the given Axes object 
    
    Parameters:
    ax : matplotlib.axes.Axes
        The Axes object on which to plot the 1D mesh.
    d : oofem.Domain
        The domain containing the mesh data.
    tstep : oofem.TimeStep
        The time step for which to plot the mesh.
    xind : int, optional
        The index of the coordinate to use for the x-axis.
    yind : int, optional
        The index of the coordinate to use for the y-axis.
    evals : list or None, optional
        The values to evaluate and color the mesh elements.
    warpField : oofem.Field, optional
        The field used to warp the mesh coordinates.
    warpScale : float, optional
        The scale factor for warping the mesh coordinates.
    label : str, optional
        The label for the colorbar.
    nodeLabels : bool, optional
        Whether to display labels for nodes.
    elementLabels : bool, optional
        Whether to display labels for elements.
    """
    x=[];
    y=[];
    
    #plt.set_title('3D OOFEM mesh')
    ax.set_aspect('equal')

    if (evals):
        # Normalize values for color mapping
        mmin = min(evals)
        mmax = max(evals)
        norm = plt.Normalize(mmin, mmax)
        colors = plt.cm.coolwarm(norm(evals))        

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
                plt.plot(x,y, color=colors[i-1], marker='s', linestyle='solid', linewidth=3, markersize=4, mec='black', mfc='black', mew=0)
            else:
                # matplotlib>3,8 will support blend_mode='screen' to solve overlapping markers issue
                plt.plot(x,y, color='black', marker='s', linestyle='solid', linewidth=0.5, markersize=3, mec='black', mfc=matplotlib.colors.to_rgba('black', 0.2), mew=0, alpha=0.3)
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
    if (evals):
        sm = plt.cm.ScalarMappable(cmap=plt.cm.coolwarm, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, label=label, ax=ax, orientation='horizontal')   
    return plt
