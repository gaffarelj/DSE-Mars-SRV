import pandas as pd
import operator
import numpy as np
import pickle

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

        zerovector = np.array([0,0,2])
        if np.dot(cp,zerovector) < 0:
            #print("Does this even work?")
            cp = -cp

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
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in new_list:
    #print(i[1],i[2],i[3],i[4][0],i[4][1],i[4][2])
    ax.quiver(i[1],i[2],i[3],i[4][0],i[4][1],i[4][2], length=0.8,arrow_length_ratio=0.6)

plt.show()