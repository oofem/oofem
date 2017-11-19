'''
Created on Aug 5, 2013

@author: carl
'''

class FEM:
    """ a simple FEM object structure"""
    def __init__(self):
        self.nnodes=0
        self.nelems=0
        self.nnodesets=0
        self.nelemsets=0
        self.nodes=[]
        self.elems=[]
        self.nodesets=[]
        self.elemsets=[]
        
class Node:
    """ a single node object """
    def __init__(self,id,coords,quadratic=1):
        self.id=id
        self.coords=coords
        self.quadratic=quadratic
        
class Element:
    """ a single finite element object"""
    def __init__(self,id,type,material,color,nnodes,cntvt):
        self.id=id
        self.type=type
        self.material=material
        self.color=color
        self.nnodes=nnodes
        self.cntvt=cntvt

class Group:
    """ FE group object """
    def __init__(self,id,name):
        self.id=id
        self.type=0
        self.name=name
        self.nitems=0
        self.items=[]

def Line2Float(line):
    """Convert a string into a list of Float"""
    return [float(x) for x in line.split()]

def Line2Int(line):
    """Convert a string into a list of Int"""
    return [int(x) for x in line.split()]
