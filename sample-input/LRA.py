import numpy as np
import matplotlib.pyplot as plt
import Nuclear as nuke
import Mesh
from Solver import Solver

#This should be the only file where things are hardcoded

mesh = Mesh.Mesh([15.0] * 11, [15.0] * 11, [30.0, 300.0, 30.0], 2)

#############################
#     Create Materials      #
#############################

F1bin = nuke.Materials(2, 1)
F1bout = nuke.Materials(2, 2)
F2bin = nuke.Materials(2, 3)
F2bout = nuke.Materials(2, 4)
Reflector = nuke.Materials(2, 5)

F1bin.Set_Chi(np.array( [1, 0] ))
F1bin.Set_D(np.array( [1.255, 0.211] ))
F1bin.Set_Sig_a(np.array( [0.008252, 0.1003] ))
F1bin.Set_NuSig_f(np.array( [0.004602, 0.1091] ))
F1bin.Set_Sig_s(np.array(( [0, 0.02533], [0, 0] )))

F1bout.Set_Chi(np.array( [1, 0] ))
F1bout.Set_D(np.array( [1.268, 0.1902] ))
F1bout.Set_Sig_a(np.array( [0.007181, 0.07047] ))
F1bout.Set_NuSig_f(np.array( [0.004609, 0.08675] ))
F1bout.Set_Sig_s(np.array(( [0, 0.02767], [0, 0] )))

F2bin.Set_Chi(np.array( [1, 0] ))
F2bin.Set_D(np.array( [1.259, 0.2091] ))
F2bin.Set_Sig_a(np.array( [0.008002, 0.08344] ))
F2bin.Set_NuSig_f(np.array( [0.004663, 0.1021] ))
F2bin.Set_Sig_s(np.array(( [0, 0.02617], [0] )))

F2bout.Set_Chi(np.array( [1, 0] ))
F2bout.Set_D(np.array( [1.259, 0.2091] ))
F2bout.Set_Sig_a(np.array( [0.008002, 0.073324] ))
F2bout.Set_NuSig_f(np.array( [0.004663, 0.1021] ))
F2bout.Set_Sig_s(np.array(( [0, 0.02617], [0, 0] )))

Reflector.Set_Chi(np.array( [1, 0] ))
Reflector.Set_D(np.array( [1.257, 0.1592] ))
Reflector.Set_Sig_a(np.array( [0.0006034, 0.01911] ))
Reflector.Set_NuSig_f(np.array( [0, 0] ))
Reflector.Set_Sig_s(np.array(( [0, 0.04754], [0, 0] )))


###############
# Set Up Mesh #
###############

q = mesh.cells

##############################
#  Assign Materials to Cells #
##############################

for j in range(len(q[0])):
	for k in range(len(q[0,j])): 
		q[0,j,k].Set_Material(Reflector)

for i in range(2):
	for j in q[1,0:2,:][i]: j.Set_Material(Reflector)
for i in range(2):
	for j in q[1,2:4,:7][i]: j.Set_Material(F2bin)
for i in range(2):
	for j in q[1,2:4,7:][i]: j.Set_Material(Reflector)
q[1,3,7].Set_Material(F2bout)
for i in q[1,4:6,0]: i.Set_Material(F1bout)
for i in range(2):
	for j in q[1,4:6,5:7][i]: j.Set_Material(F1bout)
for i in range(2):
	for j in q[1,4:6,1:5][i]: j.Set_Material(F1bin)
for i in range(2):
	for j in q[1,4:6,7:9][i]: j.Set_Material(F2bin)
for i in range(2):
	for j in q[1,4:6,9:][i]: j.Set_Material(Reflector)
for i in range(2):
	for j in q[1,6:10,:7][i]: j.Set_Material(F1bin)
for i in range(2):
	for j in q[1,8:10,:7][i]: j.Set_Material(F1bin)
for i in range(2):
	for j in q[1,6:10,7:9][i]: j.Set_Material(F2bin)
for i in range(2):
	for j in q[1,8:10,7:9][i]: j.Set_Material(F2bin)
for i in range(2):
	for j in q[1,6:10,9:][i]: j.Set_Material(Reflector)
for i in range(2):
	for j in q[1,8:10,9:][i]: j.Set_Material(Reflector)
q[1,10,0].Set_Material(F1bout)
for i in q[1,10,1:5]: i.Set_Material(F1bin)
for i in q[1,10,5:7]: i.Set_Material(F1bout)
for i in q[1,10,7:9]: i.Set_Material(F2bin)
for i in q[1,10,9:]: i.Set_Material(Reflector)

for j in range(len(q[2])):
	for k in range(len(q[2,j])): 
		q[2,j,k].Set_Material(Reflector)

###################
# Refine the Mesh #
###################

mesh = mesh.refineMesh(5.0, ['x','y'])
mesh = mesh.refineMesh(30.0, ['z'])
mesh.Set_Boundaries(["REFLECTIVE", "REFLECTIVE", "ZERO_FLUX", "ZERO_FLUX", "VACUUM", "VACUUM"])
q = mesh.cells



######################
# Plot Initial Cells #
######################

image_data = q.copy() #There are multiple layers since z was refined as well
 
for cz in range(mesh.num_cells_z):

	for cy in range(mesh.num_cells_y):
		
		for cx in range(mesh.num_cells_x):
		
			image_data[cz, cy, cx] = q[cz, cy, cx].Material.ID
			
material_data = []

for i in range(len(image_data)):
	material_data.append(image_data[i]) #Hold each z-layer in the material_data list
	
material_data = np.array(material_data, dtype='float')
	
for i in range(len(material_data)): #Name each material_data plot as the ith element of the array
	plt.imshow(material_data[i], interpolation = 'nearest', extent = [0, 165, 0, 165])
	plt.savefig("/Users/rudygarcia/Desktop/School/UROP/3D-diffusion/LRA/Material"+str(i)+".png")


##############
# Run Solver #
##############
Solver = Solver(mesh)
Solver.Compute_D()
Solver.Compute_A()
Solver.Compute_M()
Solver.Solve()

####################
# Plot Cell Fluxes #
####################

# flux0_group0 = image_data.copy() #Do this three times, once for each layer
# flux0_group1 = image_data.copy()
# flux1_group0 = image_data.copy()
# flux1_group1 = image_data.copy()
# flux2_group0 = image_data.copy()
# flux2_group1 = image_data.copy()
# 
# for cz in range(mesh.num_cells_z):
# 
# 	for cy in range(mesh.num_cells_y):
# 		
# 		for cx in range(mesh.num_cells_x):
# 		
# 			cell_id = cz * mesh.num_cells_y * mesh.num_cells_x + cy * mesh.num_cells_x + cx
# 			image_data1[cy, cx] = Solver.Phi[cell_id * 2 + 0]
# 			image_data2[cy, cx] = Solver.Phi[cell_id * 2 + 1]
# 
# image_data1 = np.array(image_data1,dtype="float")
# image_data2 = np.array(image_data2,dtype="float")
# 
# plt.imshow(image_data1,interpolation='nearest',extent=[0,165,0,165])
# 
# plt.savefig("/Users/rudygarcia/Desktop/School/UROP/3D-diffusion/LRA/Flux1.png")
# 
# plt.imshow(image_data2,interpolation='nearest',extent=[0,165,0,165])
# 
# plt.savefig("/Users/rudygarcia/Desktop/School/UROP/3D-diffusion/LRA/Flux2.png")
