""" _   __                     _         _      __        ____ ______ ____ 
   / | / /____   ____   _____ (_)____ _ (_)____/ /       /  _// ____// __ \
  /  |/ // __ \ / __ \ / ___// // __ `// // __  /______  / / / /    / /_/ /
 / /|  // /_/ // / / // /   / // /_/ // // /_/ //_____/_/ / / /___ / ____/ 
/_/ |_/ \____//_/ /_//_/   /_/ \__, //_/ \__,_/       /___/ \____//_/      
                              /____/                                             

This is an implementation of the famous N-ICP algorithm by Amberg et al. (2007).
Thanks to Yiyan Liu who provided me his Python implementation which allowed me to make this port for Rhino
and I used as a starting point for the development of my own algorithm.

This cose is a rough and messy implementation I (Giulio L.S. Sacco) made in September 2024 in an hurry. Sorry and have fun.
"""
#sys.exit(0)
#import rhinoscriptsyntax as rs
import Rhino as rn
import numpy as np
from scipy import sparse
from scipy import linalg
import math
import ghpythonlib.treehelpers as th

#debug
import sys
import time

""" DEBUG """
#GET START TIME
time_start = time.time()
time_start_string = time.strftime("%y-%m-%d %H:%M:%S", time.gmtime(time_start))

###################################################################################################
dbg = []

#PLOT Status info - Initialization
print(f"{time_start_string} UTC\n ")
with open("StatusDebug.txt", "w+") as stf:
    stf.write(f"Started at {time_start_string} UTC\n\n")

#--- --- --- ! --- --- --- SETUP --- --- --- ! --- --- --- #

#The values of alpha are assigned with an exponential decay scheme, as Yiyan Liu did and as suggested by Hasler et al., 2009
alphaValues = [maxAlpha * np.exp(-(1 / iterations * np.log(maxAlpha / minAlpha)) * i) for i in range(iterations)]
alphaValues[-1] = minAlpha #Ensure the last alpha value is exactly the minimum specified

###-### Stifness Term ###-###

# Calculate all unique edges from the source mesh triangles to aid in deformation regularisation
edges = set(tuple(sorted((face[i], face[j]))) for face in sourceMesh.Faces for i, j in [(0, 1), (0, 2), (1, 2)])
numSourceEdges = len(edges)

# M
numSourceVertexes = len(sourceMesh.Vertices)
row = []
col = []
data = []
for i, vertIndex in enumerate(edges):
    row.extend([i, i])
    col.extend([vertIndex[0], vertIndex[1]])
    data.extend([-1, 1])
matrixM = sparse.csc_matrix((data, (row, col)), shape=(numSourceEdges, numSourceVertexes))

# G exp - This is just a larger version of G. It works in the same way but it's a 12x12 bevause now X is an assembly of 12x1 instead of 4x3. 
G_exp = sparse.diags([1,1,1,gamma,1,1,1,gamma,1,1,1,gamma], shape=(12,12))

# MoxG
sparseKronMG = sparse.kron(matrixM, G_exp)

# Vector of Zeros
vecZero = np.zeros((12*numSourceEdges))

###-###  Distance Term  ###-###

# D - This D is similar to the original one, but each vertex is repeted diagonally in there lines
row = []
col = []
data = []
for i, vert in enumerate(sourceMesh.Vertices):
    for j in range(3):
        for k in range(4):
            row.append(i*3 + j)
            col.append(12*i + 4*j + k)
        data.extend([sourceMesh.Vertices[i].X, sourceMesh.Vertices[i].Y, sourceMesh.Vertices[i].Z, 1.0])
sparseD = sparse.csc_matrix((data, (row, col)), shape=(3*numSourceVertexes, 12*numSourceVertexes))

#matrixU is not here because is computed each subiterations 

###-###  Landmarks Term ###-###
if sourceLandmarks != []:
    print("Yay! You are using landmarks!\n ")

    # D_L
    numLandmark = len(sourceLandmarks)
    row = []
    col = []
    data = []
    for i, souLandInd in enumerate(sourceLandmarks):
        for j in range(3):
            for k in range(4):
                row.append(i*3 + j)
                col.append(12*souLandInd + 4*j + k)
            data.append(sourceMesh.Vertices[souLandInd].X)
            data.append(sourceMesh.Vertices[souLandInd].Y)
            data.append(sourceMesh.Vertices[souLandInd].Z)
            data.append(1.0)
    sparseDL = sparse.csc_matrix((data, (row, col)), shape=(3*numLandmark, 12*numSourceVertexes))

    #U_L
    vecU_L = []
    for vertex in targetLandmarks: 
        vecU_L.append(vertex.X)
        vecU_L.append(vertex.Y)
        vecU_L.append(vertex.Z)
    vecU_L = np.array(vecU_L)

#--- --- --- ! --- --- --- LOOP --- --- --- ! --- --- --- #

#initialize X and other stuff
X = np.tile(np.array([1,0,0,0,0,1,0,0,0,0,1,0]), numSourceVertexes)
transformedTemplate = sourceMesh.Vertices
firstLoop = True

for step, alpha in enumerate(alphaValues, 1): #OUTER LOOP
    print(f"Step {step}, alpha: {alpha}, substep:", end=" ")
    with open("StatusDebug.txt", "a") as stf:
        stf.write(f"Step {step}, alpha: {alpha}, substep:")

    # Cumpute Beta 
    if sourceLandmarks != [] and firstLoop == False:
        distanceResidual = 0
        for source, transformed in zip(transformedTemplate, sourceMesh.Vertices):
            distanceResidual += source.DistanceTo(transformed)
        
        landmarksResiduals = 0
        for i, trgLand in zip(sourceLandmarks, targetLandmarks):
            landmarksResiduals += transformedTemplate[i].DistanceTo(trgLand)

        btmp = 1.5 * (distanceResidual / landmarksResiduals)
        beta = 0 if btmp < 0 else 100 if btmp > 100 else btmp
        #print("distanceResidual", distanceResidual, "landmarksResiduals", landmarksResiduals,"Beta", beta)
 
    else: 
        firstLoop = False
        beta = 1.5        

    #Initialize the while
    oldX, subIterations = 10*X, 1
    
    while np.linalg.norm(X - oldX) >= epsilon and subIterations <= maxSubIterations: #INNER LOOP
        #print("WHILE! step", step, ",", subIterations, " - ", np.linalg.norm(X - oldX), ">=", epsilon)
        print(subIterations, end=" ")
        with open("StatusDebug.txt", "a") as stf:
            stf.write(f"{subIterations} ")

        oldX = X

        # Finding U
        target_vertices = []
        target_distances = []
        dbg = []
        if type(target) == rn.Geometry.Mesh:
            for point in transformedTemplate:
                new_target_vertices_temp = rn.Geometry.Mesh.ClosestPoint(target, point)
                dbg.append(new_target_vertices_temp)
                target_vertices.append([new_target_vertices_temp.X, new_target_vertices_temp.Y, new_target_vertices_temp.Z])
                target_distances.append(point.DistanceTo(new_target_vertices_temp))
        
        elif type(target) == rn.Geometry.PointCloud:
            for point in transformedTemplate:
                new_target_vertices_index_temp = rn.Geometry.PointCloud.ClosestPoint(target, point)
                new_target_vertices_temp = target[new_target_vertices_index_temp]
                dbg.append(new_target_vertices_temp)
                target_vertices.append([new_target_vertices_temp.X, new_target_vertices_temp.Y, new_target_vertices_temp.Z])
                target_distances.append(point.DistanceTo(new_target_vertices_temp))

        elif type(target) == rn.Geometry.Brep:
            for point in transformedTemplate:
                new_target_vertices_temp = rn.Geometry.Brep.ClosestPoint(target, point)
                dbg.append(new_target_vertices_temp)
                target_vertices.append([new_target_vertices_temp.X, new_target_vertices_temp.Y, new_target_vertices_temp.Z])
                target_distances.append(point.DistanceTo(new_target_vertices_temp))

        else:
            sys.exit("Give me proper data! I need a mesh, a proper point cloud or a brep.") # /!\ Type Hint must be on No Type Hint, otherwise it doesn't work even passing the right data!
        
        target_vertices = np.asarray(target_vertices)
        target_distances = np.asanyarray(target_distances)
        vecU = target_vertices.flatten()

        # W weights
        Wd = (target_distances <= distanceThresholdToReject).astype(int) #This is just a boolean mask, but something more refined could be theoretically by inplemented

        # W
        if Wp == []:
            W = Wd
        else:
            W = np.multiply(Wp, Wd)
        W = W.repeat(3)

        # UW vector
        vecUW = np.multiply(vecU, W)

        # W matrix multiplied W
        Wexp = sparse.kron(sparse.diags(W, shape=(3*numSourceVertexes, 3*numSourceVertexes)), np.array([1,1,1,1]), format="csc")
        sparseDW = sparseD.multiply(Wexp)

        # Assemblage of A and B matrix
        if sourceLandmarks == []:
            A = sparse.vstack([alpha*sparseKronMG, sparseDW])
            B = np.concatenate((vecZero, vecUW))
        else:
            A = sparse.vstack([alpha*sparseKronMG, sparseDW, beta*sparseDL]) #In equation 12 (Amberg et al.2007) beta is only on D_L but according to equation 5 i assumes it was intended to be also on U_L.  
            B = np.concatenate((vecZero, vecUW, beta*vecU_L))

        # Solve
        X, istop, itn, normt = sparse.linalg.lsqr(A, B, atol=0, btol=0, conlim=0, iter_lim=10000)[:4]
        if itn > 9999:
            print("Ok, i've been lazy here and I skipped the study of sparse.linalg.lsqr leveraging the power of my pc and my usecase. If you are reading this message you should probably take a look here.")
        #print("Solution", istop, itn, normt)

        # Transform template (some post projects at subiteration level)
        transformedTemplateTemp = sparseD @ X
        transformedTemplateTemp = transformedTemplateTemp.reshape(-1, 3)
        transformedTemplate = []
        for point in transformedTemplateTemp:
            transformedTemplate.append(rn.Geometry.Point3d(point[0], point[1], point[2]))

        subIterations += 1
    
    print("")
    with open("StatusDebug.txt", "a") as stf:
        stf.write(f"\n")

#--- --- --- ! --- --- --- POST PROCESSING --- --- --- ! --- --- --- #

#Creation of the registeredMesh in Rhino format
registeredMesh = rn.Geometry.Mesh()#Initialization
for vrt in transformedTemplate: 
    registeredMesh.Vertices.Add(vrt)
for fc in sourceMesh.Faces:   
    registeredMesh.Faces.AddFace(fc)

# displacement field in Rhino format
displacementField = []
for start, finish in zip(sourceMesh.Vertices, registeredMesh.Vertices):
    displacementField.append(rn.Geometry.Vector3d(finish - start))
displacementField = th.list_to_tree(displacementField, source=[0,0]) #this th.list_to_tree is black magic to me.

#--- --- --- ! --- --- --- END --- --- --- ! --- --- --- #

#FINAL TIME TAG
time_end = time.time()
time_end_string = time.strftime("%y-%m-%d %H:%M:%S", time.gmtime(time_end))
time_elapsed = time_end - time_start 
time_elapsed_minutes = time_elapsed // 60 
time_elapsed_seconds = round(time_elapsed % 60, 2)
print(f" \nDONE AT {time_end_string} UTC!!!\n \nElapsed time: {time_elapsed_minutes} minutes and {time_elapsed_seconds} seconds")
with open("StatusDebug.txt", "a") as stf:
    stf.write(f"\nDONE AT {time_end_string} UTC!!!\n \nElapsed time: {time_elapsed_minutes} minutes and {time_elapsed_seconds} seconds") 

print(" \n \n \n \n \n \n \n \n \nฅ^•ﻌ•^ฅ")
