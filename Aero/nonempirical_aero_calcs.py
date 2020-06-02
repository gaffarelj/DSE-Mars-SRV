import pandas as pd
import operator
import numpy as np
import matplotlib
import matplotlib.cm as cmx
import pickle
import astrodynamics.mars_standard_atmosphere as MSA



mesh_elements = pd.read_csv("surfacemesh.csv",sep=";")

mesh_elements = pd.DataFrame(mesh_elements).to_numpy()

def mesh_elements_with_vector(mesh_elements):
    uber_uber_list = []
    for i in range(len(mesh_elements)):
        xloc = mesh_elements[i][1]
        yloc = mesh_elements[i][2]
        zloc = mesh_elements[i][3]

        distancedict = {}
        for j in range(len(mesh_elements)):
            elementnumber = mesh_elements[j][0]
            xloc2 = mesh_elements[j][1]
            yloc2 = mesh_elements[j][2]
            zloc2 = mesh_elements[j][3]
            distance = ((xloc-xloc2)**2+(yloc-yloc2)**2+(zloc-zloc2)**2)**0.5

            distancedict[elementnumber] = distance

        sorted_tuple = sorted(distancedict.items(), key = operator.itemgetter(1))


        uber_list = [mesh_elements[i][0],xloc,yloc,zloc]
        vector_list = []
        for i in range(3):
            node_number = list(sorted_tuple[i+1])[0]
            vector_list.append((mesh_elements[int(node_number)-1][-3:]))

        p1 = vector_list[0]
        p2 = vector_list[1]
        p3 = vector_list[2]

        v1 = p3-p1
        v2 = p2 - p1
        cp = np.cross(v1,v2)

        zerovector = np.array([xloc,yloc,zloc])
        if np.dot(cp,zerovector) > 0:
            #print("Does this even work?")
            cp = -cp

        cp = cp / np.linalg.norm(cp)


        uber_list.append(cp)
        uber_uber_list.append(uber_list)

    with open("mesh_with_vectors.txt", "wb") as f:
        pickle.dump(uber_uber_list,f)
    return uber_uber_list

#some_list = mesh_elements_with_vector(mesh_elements)

with open("mesh_with_vectors.txt", "rb") as f:
  new_list = pickle.load(f)

#print(new_list)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

#for i in new_list:
    #print(i[1],i[2],i[3],i[4][0],i[4][1],i[4][2])
#    ax.quiver(i[1],i[2],i[3],i[4][0],i[4][1],i[4][2], length=0.8,arrow_length_ratio=0.6)





def aero_angles(alpha, mesh_list):

    '''
    Takes in alpha in deg, and mesh_list as it is set in mesh_elements_with_vectors().
    '''
    alpha = alpha*np.pi/180

    #print(np.cos(alpha),np.sin(alpha))

    aero_vector = np.array([0, np.sin(alpha),-np.cos(alpha)])
    uptream_mesh_list = []
    downstream_mesh_list = []
    for i in range(len(mesh_list)):
        mesh_vector = mesh_list[i][4]
        #mesh_vector[0] = 0

        dotprod = np.dot(aero_vector,mesh_vector)
        mesh_item = mesh_list[i]
        if dotprod > 0.05:
            angle1 = np.arctan(aero_vector[1]/aero_vector[2])
            angle2 = np.arctan(mesh_vector[1]/mesh_vector[2])
            angle = angle1-angle2
            mesh_item.append(angle)
            uptream_mesh_list.append(mesh_item)
        if dotprod <= 0.05:
            downstream_mesh_list.append(mesh_item)

    #for i in uptream_mesh_list:
        #print(i[1],i[2],i[3],i[4][0],i[4][1],i[4][2])
    #    ax.quiver(i[1],i[2],i[3],i[4][0],i[4][1],i[4][2], length=0.8,arrow_length_ratio=0.6)
        #ax.scatter(i[1],i[2],i[3])
    #ax.quiver(0, 0, 16, 0, np.sin(alpha), -np.cos(alpha), color='r')
    return uptream_mesh_list, downstream_mesh_list

def aero_cp_calcs(upstream_mesh,downstream_mesh, p_inf, p_02,M_inf,alpha):
    gamma = 1.37
    Cp_max = 2/(gamma*M_inf**2)*(p_02/p_inf -1)
    #print(Cp_max)
    perp_alpha = (alpha + 90) * np.pi/180
    alpha = alpha * np.pi/180
    aero_vector = np.array([0, np.sin(alpha), -np.cos(alpha)])
    perp_aero_vector = np.array([0, np.sin(perp_alpha), -np.cos(perp_alpha)])

    pressure_list = []
    areasum = 0
    for i in range(len(upstream_mesh)):
        mesh_vector = upstream_mesh[i][4]
        mesh_vector = mesh_vector/np.linalg.norm(mesh_vector)
        area = abs(np.dot(aero_vector,mesh_vector))/(len(upstream_mesh))
        area = area * (len(upstream_mesh)+len(downstream_mesh))/len(upstream_mesh)
        angle = upstream_mesh[i][5]
        Cp = Cp_max*(np.cos(angle)**2)
        Cp = Cp*area
        areasum = areasum + area
        mesh_vector = mesh_vector * Cp
        pressure_list_item = [upstream_mesh[i][1],upstream_mesh[i][2],upstream_mesh[i][3],Cp, mesh_vector]
        pressure_list.append(pressure_list_item)
        #print(Cp)
    for i in range(len(downstream_mesh)):
        mesh_vector = np.array([0,0,0])
        pressure_list_item = [downstream_mesh[i][1],downstream_mesh[i][2],downstream_mesh[i][3],0, mesh_vector]
        pressure_list.append(pressure_list_item)
    #print("AREA", areasum)
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')

    vectorsum = np.array([0,0,0])
    for i in pressure_list:
        vectorsum = vectorsum + i[4]

        i[4][0] = i[4][0]
        i[4][1] = i[4][1]
        i[4][2] = i[4][2]
        #ax.quiver(i[0], i[1], i[2], i[4][0], i[4][1], i[4][2], arrow_length_ratio=0.6)

    vectorsum = vectorsum
    aero_vector = aero_vector
    #print("Cp:")
    #print(np.linalg.norm(vectorsum))

    perp_aero_vector = perp_aero_vector

    projected_vector = aero_vector * np.dot(vectorsum,aero_vector)/np.dot(aero_vector,aero_vector)
    projected_vector_perp = perp_aero_vector * np.dot(vectorsum,perp_aero_vector)/np.dot(perp_aero_vector,perp_aero_vector)
    #print("Cd,Cl:")
    cd = np.linalg.norm(projected_vector)
    cl = np.linalg.norm(projected_vector_perp)
    #print((cl**2 + cd **2)**0.5)
    #print(np.linalg.norm(projected_vector),np.linalg.norm(projected_vector_perp))
    #print(projected_vector,projected_vector_perp)
    xlist = [16,16,-16,-16]
    ylist = [16,-16,-16,16]
    #ax.scatter(xlist,ylist)

    #ax.quiver(0,0,0,aero_vector[0],aero_vector[1],aero_vector[2], color = 'g')
    #ax.quiver(0,0,7,vectorsum[0],vectorsum[1],vectorsum[2])
    #ax.quiver(0,0,7,projected_vector[0],projected_vector[1],projected_vector[2],color = 'r')
    #ax.quiver(0,0,7,projected_vector_perp[0],projected_vector_perp[1],projected_vector_perp[2],color = 'r')
    #ax.quiver(0,0,0,perp_aero_vector[0],perp_aero_vector[1],perp_aero_vector[2],color = 'g')

    #plt.show()

    #print(vectorsum)


    return cl, cd

def ultimate_aero_functions(p_inf, p_02, M_inf, alpha):
    with open("mesh_with_vectors.txt", "rb") as f:
        new_list = pickle.load(f)

    upstream_mesh, downstream_mesh = aero_angles(alpha, new_list)

    cl, cd = aero_cp_calcs(upstream_mesh, downstream_mesh, p_inf, p_02, M_inf, alpha)

    return cl, cd

p_inf = MSA.get_pressure(20000)
p_02 = 783
M_inf = 3
alpha = 50

#upstream_mesh, downstream_mesh = aero_angles(alpha,new_list)

#pressures = aero_cp_calcs(upstream_mesh,downstream_mesh,p_inf,p_02,M_inf,alpha)




cl, cd = ultimate_aero_functions(p_inf,p_02, M_inf, alpha)
print(cl, cd)


xlist = []
ylist = []
zlist = []
cplist = []
#for i in pressures:
#    xlist.append(i[0])
#    ylist.append(i[1])
#    zlist.append(i[2])
#    cplist.append(i[3])

xlist = np.array(xlist)
ylist = np.array(ylist)
zlist = np.array(zlist)
cplist = np.array(cplist)



def scatter3d(x,y,z, cs, colorsMap='jet'):
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, c=scalarMap.to_rgba(cs))
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    plt.show()

#scatter3d(xlist,ylist,zlist,cplist)