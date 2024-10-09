lx = 1.
ly = 1.
lz = 1.
nx = 10
ny = 10
nz = 10
dx = lx/nx
dy = ly/ny
dz = lz/nz


nodes = []
elements = []
nn = 0
ne = 0

xsurfoffset = 0
# generate external boundary entities first
# x_surfs
for k in range(nz):
    for j in range(ny):
        y=j*dy
        z=k*dz
        # front face
        nodes.append((nn+1,0, y,z))
        nodes.append((nn+2,0, y+dy,z))
        nodes.append((nn+3,0, y+dy,z+dz))
        nodes.append((nn+4,0, y,z+dz))
        nn+=4
        # rear face
        nodes.append((nn+1,lx, y,z))
        nodes.append((nn+2,lx, y+dy,z))
        nodes.append((nn+3,lx, y+dy,z+dz))
        nodes.append((nn+4,lx, y,z+dz))
        nn+=4
# y surfs
ysurfoffset = nn
for k in range(nz):
    for i in range(nx):
        x=i*dx
        z=k*dz
        # front face
        nodes.append((nn+1,x, 0,z+dz))
        nodes.append((nn+2,x+dx, 0,z+dz))
        nodes.append((nn+3,x+dx, 0,z))
        nodes.append((nn+4,x, 0,z))
        nn+=4
        # rear face
        nodes.append((nn+1,x, ly,z+dz))
        nodes.append((nn+2,x+dx, ly,z+dz))
        nodes.append((nn+3,x+dx, ly,z))
        nodes.append((nn+4,x, ly,z))
        nn+=4
# z surfs
zsurfoffset = nn
for j in range(ny):
    for i in range(nx):
        x=i*dx
        y=j*dy
        # front face
        nodes.append((nn+1,x, y,0))
        nodes.append((nn+2,x+dx, y,0))
        nodes.append((nn+3,x+dx, y+dy,0))
        nodes.append((nn+4,x, y+dy,0))
        nn+=4
        # rear face
        nodes.append((nn+1,x, y,lz))
        nodes.append((nn+2,x+dx, y,lz))
        nodes.append((nn+3,x+dx, y+dy,lz))
        nodes.append((nn+4,x, y+dy,lz))
        nn+=4

# generate internal nodes & elements      

for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            # generate segment nodes
            x = i*dx
            y = j*dy
            z = k*dz
            nodes.append((nn+1,x,y,z))
            nodes.append((nn+2,x+dx,y,z))
            nodes.append((nn+3,x+dx,y+dy,z))
            nodes.append((nn+4,x,y+dy,z))
            nodes.append((nn+5,x,y,z+dz))
            nodes.append((nn+6,x+dx,y,z+dz))
            nodes.append((nn+7,x+dx,y+dy,z+dz))
            nodes.append((nn+8,x,y+dy,z+dz))
            
            # generate element
            elements.append(("sadgbrick1", ne+1, 8, nn+1, nn+2, nn+3, nn+4, nn+5, nn+6, nn+7, nn+8))
            ne+=1
            # generate interfaces
            if i == 0:
                indx = xsurfoffset+((k)*ny+j)*8
                elements.append(("sadgbquad1", ne+1, 8, indx+1, indx+2, indx+3, indx+4, nn+1, nn+4, nn+8, nn+5))
                ne+=1
            else:
                indx = nn-8
                elements.append(("sadgBquad1", ne+1, 8, indx+2, indx+3, indx+7, indx+6, nn+1, nn+4, nn+8, nn+5))
                ne+=1

            if i == nx-1:
                indx = xsurfoffset+((k)*ny+j)*8+4
                elements.append(("sadgbquad1", ne+1, 8, nn+2, nn+3, nn+7, nn+6, indx+1, indx+2, indx+3, indx+4))
                ne+=1
            
            if j == 0:
                indx = ysurfoffset+((k)*nx+i)*8
                elements.append(("sadgbquad1", ne+1, 8, indx+1, indx+2, indx+3, indx+4, nn+5, nn+6, nn+2, nn+1))
                ne+=1
            else:
                indx = nn-nx*8
                elements.append(("sadgBquad1", ne+1, 8, indx+8, indx+7, indx+3, indx+4, nn+5, nn+6, nn+2, nn+1))
                ne+=1
            if j == ny-1:
                indx = ysurfoffset+((k)*nx+i)*8+4
                elements.append(("sadgbquad1", ne+1, 8, nn+8, nn+7, nn+3, nn+4, indx+1, indx+2, indx+3, indx+4))
                ne+=1

            if k == 0:
                indx = zsurfoffset+((j)*nx+i)*8
                elements.append(("sadgbquad1", ne+1, 8, indx+1, indx+2, indx+3, indx+4, nn+1, nn+2, nn+3, nn+4))
                ne+=1 
            else:
                indx = nn-nx*ny*8
                elements.append(("sadgBquad1", ne+1, 8, indx+5, indx+6, indx+7, indx+8, nn+1, nn+2, nn+3, nn+4))
                ne+=1
            if k == nz-1:
                indx = zsurfoffset+((j)*nx+i)*8+4
                elements.append(("sadgbquad1", ne+1, 8, nn+5, nn+6, nn+7, nn+8, indx+1, indx+2, indx+3, indx+4))
                ne+=1
            nn+=8

print("nnode %d nelem %d"%(len(nodes),len(elements)))
#print nodes
for n in nodes:
    print ("Node %5d coords 3 %12.5e %12.5e %12.5e"%(n[0],n[1],n[2],n[3]))
#print elements
for e in elements:
    print ("%10s %4d nodes %2d %5d %5d %5d %5d %5d %5d %5d %5d"%(e[0],e[1],e[2],e[3],e[4],e[5],e[6],e[7],e[8],e[9],e[10]))
# external boundary set\
print ("set 2 noderanges {(1 %d)}"%(8*(nx*ny+nx*nz+ny*nz),))
