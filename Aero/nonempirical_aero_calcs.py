import pandas as pd
import operator
import numpy as np
import matplotlib
import matplotlib.cm as cmx
import pickle
import astrodynamics.mars_standard_atmosphere as MSA
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def mesh_elements_with_vector(input_filename,output_filename,gamma):
    mesh_elements = pd.read_csv(input_filename, sep=";")
    gamma = gamma*np.pi/180
    mesh_elements = pd.DataFrame(mesh_elements).to_numpy()

    for i in range(len(mesh_elements)):
        mesh_elements[i][2] = mesh_elements[i][2] + np.sin(gamma) * mesh_elements[i][3]
        mesh_elements[i][3] = mesh_elements[i][3] - np.sin(gamma) * mesh_elements[i][3]

    uber_uber_list = []
    for i in range(len(mesh_elements)):

        xloc = mesh_elements[i][1]
        yloc = mesh_elements[i][2]#+np.sin(gamma)*mesh_elements[i][3]
        zloc = mesh_elements[i][3]#-0.93#-np.sin(gamma)*mesh_elements[i][3]

        distancedict = {}
        for j in range(len(mesh_elements)):

            elementnumber = j #mesh_elements[j][0]
            xloc2 = mesh_elements[j][1]
            yloc2 = mesh_elements[j][2]#+np.sin(gamma)*mesh_elements[j][3]
            zloc2 = mesh_elements[j][3]#-0.93#-np.sin(gamma)*mesh_elements[j][3]
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
            cp = -cp

        cp = cp / np.linalg.norm(cp)

        uber_list.append(cp)
        uber_uber_list.append(uber_list)
        #print(cp)

    with open(output_filename, "wb") as f:
        pickle.dump(uber_uber_list,f)

    #print("Succesfully converted mesh to txt!")
    return

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
        if dotprod > 0.05 and mesh_vector[2] is not 0:
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
        #print(mesh_vector)
        mesh_vector = mesh_vector * Cp

        if upstream_mesh[i][2] < 5.5:
            pressure_list_item = [upstream_mesh[i][1],upstream_mesh[i][2],upstream_mesh[i][3],Cp, mesh_vector]
        if upstream_mesh[i][2] >= 5.5:
            pressure_list_item = [upstream_mesh[i][1], upstream_mesh[i][2], upstream_mesh[i][3], 0, mesh_vector]
        pressure_list.append(pressure_list_item)
        #print(pressure_list_item)
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

    cplist = []
    plist = []

    for i in pressure_list:
        #print(i)
        cplist.append(i[2]*i[3])
        plist.append(i[3])

    #print(sum(plist))
    cp = sum(cplist)/sum(plist)
    print("Cp:", cp)

    #ax.quiver(0,0,cp,aero_vector[0],aero_vector[1],aero_vector[2], color = 'g')
    #ax.quiver(0,0,cp,vectorsum[0],vectorsum[1],vectorsum[2])
    #ax.quiver(0,0,7,projected_vector[0],projected_vector[1],projected_vector[2],color = 'r')
    #ax.quiver(0,0,7,projected_vector_perp[0],projected_vector_perp[1],projected_vector_perp[2],color = 'r')
    #ax.quiver(0,0,cp,perp_aero_vector[0],perp_aero_vector[1],perp_aero_vector[2],color = 'g')

    #plt.show()


    #print(vectorsum)

    cgpos = 9.24
    clmoment = cl*np.cos(alpha)*(cgpos-cp)
    cdmoment = cd*np.sin(alpha)*(cgpos-cp)

    rho = 0.00555

    aeromoment = (clmoment+cdmoment)*0.5*rho*553**2 * 350

    #print("Aerodynamic moment:", aeromoment)

    return cl, cd, pressure_list


def ultimate_aero_functions(filename, p_inf, p_02, M_inf, alpha):
    with open(filename, "rb") as f:
        new_list = pickle.load(f)

    upstream_mesh, downstream_mesh = aero_angles(alpha, new_list)

    cl, cd, pressures = aero_cp_calcs(upstream_mesh, downstream_mesh, p_inf, p_02, M_inf, alpha)

    return cl, cd, pressures


def pitching_control_functions(filename,filename_BF, p_inf, p_02, M_inf, alpha,gamma):
    with open(filename, "rb") as f:
        new_list = pickle.load(f)
    with open(filename_BF, "rb") as f:
        BF_list = pickle.load(f)

    for i in range(len(BF_list)):
        new_list.append(BF_list[i])

    upstream_mesh, downstream_mesh = aero_angles(alpha, new_list)

    cl, cd, pressures = aero_cp_calcs(upstream_mesh, downstream_mesh, p_inf, p_02, M_inf, alpha)

    return cl, cd, pressures


def aero_normalshock(altitude, mach_initial):
    altitude = altitude
    gamma = MSA.gamma
    t_static = MSA.get_temperature(altitude)
    p_static = MSA.get_pressure(altitude)
    rho_static = MSA.get_density(p_static, t_static)

    mach = np.sqrt((1 + ((gamma - 1) / 2 * mach_initial ** 2)) / (gamma * mach_initial ** 2 - (gamma - 1) / 2))
    pressure = p_static * (1 + 2 * gamma / (gamma + 1) * (mach_initial ** 2 - 1))
    density = rho_static * (((gamma + 1) * mach_initial ** 2) / (2 + (gamma - 1) * mach_initial ** 2))
    temperature = t_static * pressure / p_static * rho_static / density
    a = np.sqrt(gamma * MSA.R * temperature)

    V = mach * a
    dyn_pressure = 0.5 * density * V * V
    return pressure, dyn_pressure

def write_aero_to_csv(filename):
    alphalist = np.arange(30,60.1,0.1)
    machlist = np.arange(0.1,20.1,0.1)
    cl_matrix = np.zeros((len(machlist)+1,len(alphalist)+1))
    cd_matrix = np.zeros((len(machlist)+1,len(alphalist)+1))
    cl_matrix[0][0] = 0
    cd_matrix[0][0] = 0
    for j in range(len(machlist)):
        cl_matrix[j+1][0] = machlist[j]
        cd_matrix[j+1][0] = machlist[j]

        for i in range(len(alphalist)):
            cl_matrix[0][i+1] = alphalist[i]
            cd_matrix[0][i+1] = alphalist[i]
            p_inf = MSA.get_pressure(10000)
            p_02, irrelevant = aero_normalshock(10000, machlist[j])
            cl,cd,pressures = ultimate_aero_functions(filename,p_inf,p_02,machlist[j],alphalist[i])
            cl_matrix[j+1][i+1] = cl
            cd_matrix[j+1][i+1] = cd

        print("Progress:",(j+1)/len(machlist)*100,"%")

    with open ('cl_standard_config.csv', mode = 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter = ';', quotechar= '"', quoting=csv.QUOTE_MINIMAL)
        #csv_writer.writerow(alphalist)
        for i in range(len(cl_matrix)):
            csv_writer.writerow(cl_matrix[i])

    with open ('cd_standard_config.csv', mode = 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter = ';', quotechar= '"', quoting=csv.QUOTE_MINIMAL)
        #csv_writer.writerow(alphalist)
        for i in range(len(cd_matrix)):
            csv_writer.writerow(cd_matrix[i])

    return


'''
alphalist = np.arange(10,70,10)
machlist = np.arange(5,25,5)
cl_matrix = np.zeros((len(machlist)+1,len(alphalist)+1))
cd_matrix = np.zeros((len(machlist)+1,len(alphalist)+1))
cl_matrix[0][0] = 0
cd_matrix[0][0] = 0

for j in range(len(machlist)):
    cl_matrix[j+1][0] = machlist[j]
    cd_matrix[j+1][0] = machlist[j]

    for i in range(len(alphalist)):
        cl_matrix[0][i+1] = alphalist[i]
        cd_matrix[0][i+1] = alphalist[i]
        p_inf = MSA.get_pressure(10000)
        p_02, irrelevant = aero_normalshock(10000, machlist[j])
        cl,cd,pressures = ultimate_aero_functions(p_inf,p_02,machlist[j],alphalist[i])
        cl_matrix[j+1][i+1] = cl
        cd_matrix[j+1][i+1] = cd

    print("Progress:",(j+1)/len(machlist)*100,"%")
print(cd_matrix)'''

#upstream_mesh, downstream_mesh = aero_angles(alpha,new_list)

#pressures = aero_cp_calcs(upstream_mesh,downstream_mesh,p_inf,p_02,M_inf,alpha)

#Vel = 533m/s
#Altitude = 8.6km
#



M_inf = 3
alpha = 50
p_inf = MSA.get_pressure(8600)
p_02, irrelevant = aero_normalshock(8600,M_inf)


p_02 = 1642.8
p_inf = 103.86
M_inf = 3.316
vel = 553.2

#74585.8 Nm

#mesh_elements_with_vector("009_surface_mesh.csv","009_aeroshell_mesh.txt",0)

gamma= 10
gammalist = [10]
for i in range(len(gammalist)):
    gamma = gammalist[i]
    print("Gamma:",gamma)
#mesh_elements_with_vector("008_surface_mesh.csv","008_aeroshell_mesh.txt",0)
    mesh_elements_with_vector("body_flap_mesh2.csv","body_flap_mesh.txt",gamma)
    cl, cd, pressures = pitching_control_functions("009_aeroshell_mesh.txt","body_flap_mesh.txt",p_inf,p_02, M_inf, alpha,gamma)

#print(cl, cd)

cplist = []
plist = []

for i in pressures:
    #print(i)
    cplist.append(i[2])
    plist.append(i[3])

cp = sum(cplist)/len(pressures)
#print("Cp z-loc: ", cp)

#print(len(pressures))


xlist = []
ylist = []
zlist = []
cplist = []
for i in pressures:
    if i[3]>0.0005:
        xlist.append(i[0])
        ylist.append(i[1])
        zlist.append(i[2])
        cplist.append(i[3])

xlist = np.array(xlist)
ylist = np.array(ylist)
zlist = np.array(zlist)
cplist = np.array(cplist)

def scatter3d(x,y,z, cs, colorsMap='jet'):
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=0, vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    #z = np.array([z.T,z])
    #ax.plot_wireframe(x,y,z)
    #ax.scatter(x, y, z, c=scalarMap.to_rgba(cs))

    xlist2, ylist2 = np.meshgrid(xlist, ylist)
    zlist2 , zlist3 = np.meshgrid(zlist, zlist)
    print(xlist2)
    print(ylist2)
    #print(zlist2)
    #print(zlist3)
    ax.plot_surface(xlist2,ylist2,zlist2,rstride=1,cstride= 1, cmap = cm)
    #ax.scatter(xlist2,ylist2)
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    plt.show()

#scatter3d(xlist,ylist,zlist,cplist)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)

print(X)
print(Y)

R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)

print(Z)

fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.viridis)

plt.show()