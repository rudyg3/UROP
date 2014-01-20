import numpy as np
import matplotlib.pyplot as plt
import Nuclear as nuke
import Mesh
from math import *
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import *
import time

class Solver(object):
	def __init__(self, mesh): #mesh being an instance of the Mesh class
		self.mesh = mesh
		self.num_groups = self.mesh.num_groups #mesh.cells[0,0,0].Material.num_groups
		self.A = lil_matrix((mesh.num_cells * self.num_groups , mesh.num_cells * self.num_groups))
		self.M = lil_matrix((mesh.num_cells * self.num_groups , mesh.num_cells * self.num_groups))
		self.Phi = np.ones(mesh.num_cells*self.num_groups)
		self.keff = 1.0
		self.keff_old = 1.0
	
	def Compute_D(self):
		print "Computing D values for surfaces..."
		self.d_start = time.time()
		lengths_x = self.mesh.lengths_x
		lengths_y = self.mesh.lengths_y
		lengths_z = self.mesh.lengths_z
		
		for z in range(self.mesh.num_cells_z):
			for y in range(self.mesh.num_cells_y):
				for x in range(self.mesh.num_cells_x):
					cell = self.mesh.cells[z,y,x]

					for face in range(6):
						for group in range(self.num_groups):
							if face == 0 and x != 0:
								cell_next = self.mesh.cells[z,y,x-1]
								D = cell.Material.D[group]
								D_next = cell_next.Material.D[group]
								D_surface = ( 2 * D * D_next ) / (D * lengths_x[x] + D_next * lengths_x[x-1])
								cell.D_surface[group,0] = D_surface
							
							elif face == 0 and x == 0:
								D = cell.Material.D[group]
								
								if self.mesh.boundary[0] == "VACUUM":
									D_surface = (2 * D) / lengths_x[x] / (1 + 4 * D / lengths_x[x]) 
									print D_surface
									cell.D_surface[group,0] = D_surface
									
								
								elif self.mesh.boundary[0] == "ZERO_FLUX":
									D_surface = (2 * D) / lengths_x[x]
									cell.D_surface[group,0] = D_surface
								
								elif self.mesh.boundary[0] == "REFLECTIVE":
									cell.D_surface[group,0] = 0.0
								
	
							elif face == 1 and (y != self.mesh.num_cells_y-1):
								cell_next = self.mesh.cells[z,y+1,x]
								D = cell.Material.D[group]
								D_next = cell_next.Material.D[group]
								D_surface = ( 2 * D * D_next ) / (D * lengths_y[y] + D_next * lengths_y[y+1])
								cell.D_surface[group,1] = D_surface
								
							elif face == 1 and not (y != self.mesh.num_cells_y-1):
								D = cell.Material.D[group]
								
								if self.mesh.boundary[1] == "VACUUM":
									D_surface = (2 * D) / lengths_y[y] / (1 + 4 * D / lengths_y[y]) 
									cell.D_surface[group,1] = D_surface
								
								elif self.mesh.boundary[1] == "ZERO_FLUX":
									D_surface = (2 * D) / lengths_y[y]
									cell.D_surface[group,1] = D_surface
								
								elif self.mesh.boundary[1] == "REFLECTIVE":
									cell.D_surface[group,1] = 0.0

							elif face == 2 and (x != self.mesh.num_cells_x-1):
								cell_next = self.mesh.cells[z,y,x+1]
								D = cell.Material.D[group]
								D_next = cell_next.Material.D[group]
								D_surface = ( 2 * D * D_next ) / (D * lengths_x[x] + D_next * lengths_x[x+1])
								cell.D_surface[group,2] = D_surface
								
							elif face == 2 and not (x != self.mesh.num_cells_x-1):
								D = cell.Material.D[group]
								
								if self.mesh.boundary[2] == "VACUUM":
									D_surface = (2 * D) / lengths_x[x] / (1 + 4 * D / lengths_x[x]) 
									cell.D_surface[group,2] = D_surface
								
								elif self.mesh.boundary[2] == "ZERO_FLUX":
									D_surface = (2 * D) / lengths_x[x]
									cell.D_surface[group,2] = D_surface
								
								elif self.mesh.boundary[2] == "REFLECTIVE":
									cell.D_surface[group,2] = 0.0
	
							elif face == 3 and (y != 0):
								cell_next = self.mesh.cells[z,y-1,x]
								D = cell.Material.D[group]
								D_next = cell_next.Material.D[group]
								D_surface = ( 2 * D * D_next ) / (D * lengths_y[y] + D_next * lengths_y[y-1])
								cell.D_surface[group,3] = D_surface
								
							elif face == 3 and not (y != 0):
								D = cell.Material.D[group]
								
								if self.mesh.boundary[3] == "VACUUM":
									D_surface = (2 * D) / lengths_y[y] / (1 + 4 * D / lengths_y[y]) 
									cell.D_surface[group,3] = D_surface
								
								elif self.mesh.boundary[3] == "ZERO_FLUX":
									D_surface = (2 * D) / lengths_y[y]
									cell.D_surface[group,3] = D_surface
								
								elif self.mesh.boundary[3] == "REFLECTIVE":
									cell.D_surface[group,3] = 0.0
								
							elif face == 4 and (z != 0):
								cell_next = self.mesh.cells[z-1,y,x]
								D = cell.Material.D[group]
								D_next = cell_next.Material.D[group]
								D_surface = ( 2 * D * D_next ) / (D * lengths_z[z] + D_next * lengths_z[z-1])
								cell.D_surface[group,4] = D_surface
								
							elif face == 4 and not (z != 0):
								D = cell.Material.D[group]
								
								if self.mesh.boundary[4] == "VACUUM":
									D_surface = (2 * D) / lengths_z[z] / (1 + 4 * D / lengths_z[z]) 
									cell.D_surface[group,4] = D_surface
								
								elif self.mesh.boundary[4] == "ZERO_FLUX":
									D_surface = (2 * D) / lengths_z[z]
									cell.D_surface[group,4] = D_surface
								
								elif self.mesh.boundary[4] == "REFLECTIVE":
									cell.D_surface[group,4] = 0.0
							
							elif face == 5 and (z != self.mesh.num_cells_z - 1):
								cell_next = self.mesh.cells[z+1,y,x]
								D = cell.Material.D[group]
								D_next = cell_next.Material.D[group]
								D_surface = ( 2 * D * D_next ) / (D * lengths_z[z] + D_next * lengths_z[z+1])
								cell.D_surface[group,5] = D_surface
								
							elif face == 5 and not (z != self.mesh.num_cells_z - 1):
								D = cell.Material.D[group]
								
								if self.mesh.boundary[5] == "VACUUM":
									D_surface = (2 * D) / lengths_z[z] / (1 + 4 * D / lengths_z[z]) 
									cell.D_surface[group,5] = D_surface
								
								elif self.mesh.boundary[5] == "ZERO_FLUX":
									D_surface = (2 * D) / lengths_z[z]
									cell.D_surface[group,5] = D_surface
								
								elif self.mesh.boundary[5] == "REFLECTIVE":
									cell.D_surface[group,5] = 0.0
									
									
		self.d_end = time.time()
	

###############################
# Get Initial Matrix Equation #
###############################
	
	def Compute_A(self):
		print "Creating loss matrix..."
		self.a_start = time.time()
		ng = self.num_groups
		cx = self.mesh.num_cells_x
		cy = self.mesh.num_cells_y
		cz = self.mesh.num_cells_z
		lengths_x = self.mesh.lengths_x
		lengths_y = self.mesh.lengths_y
		lengths_z = self.mesh.lengths_z
		
		self.A = self.A.tolil()
		self.A = self.A * 0.0
		
		for z in range(cz):
			for y in range(cy):
				for x in range(cx):
					cell = self.mesh.cells[z,y,x]
					cell_id = (z * cx * cy + y * cx + x)

					for g in range(ng):
					
						#Removal
						self.A[cell_id * ng + g, cell_id * ng + g] = cell.volume * cell.Material.Sig_a[g]
						
						for e in range(ng):
							if e != g:
								self.A[cell_id * ng + g, cell_id * ng + e] = -1 * cell.volume * cell.Material.Sig_s[e][g]
								self.A[cell_id * ng + g, cell_id * ng + g] += cell.volume * cell.Material.Sig_s[g][e]

						#Transport to left cell
						if x != 0:
							self.A[cell_id * ng + g, (cell_id - 1) * ng + g] = -1 * lengths_y[y] * lengths_z[z] * cell.D_surface[g][0]
						
						self.A[cell_id * ng + g, cell_id * ng + g] += lengths_y[y] * lengths_z[z] * cell.D_surface[g][0]
					
				
						#Transport to south cell
						if y != cy - 1:
							self.A[cell_id * ng + g, (cell_id + cx) * ng + g] = -1 * lengths_x[x] * lengths_z[z] * cell.D_surface[g][1]
							
						self.A[cell_id * ng + g, cell_id * ng + g] += lengths_x[x] * lengths_z[z] * cell.D_surface[g][1]

						#Transport to right cell
						if x != cx-1:
							self.A[cell_id * ng + g, (cell_id + 1) * ng + g] = -1 * lengths_y[y] * lengths_z[z] * cell.D_surface[g][2]
						
						self.A[cell_id * ng + g, cell_id * ng + g] += lengths_y[y] * lengths_z[z] * cell.D_surface[g][2]
	
						#Transport to north cell
						if y != 0:
							self.A[cell_id * ng + g, (cell_id - cx) * ng + g] = -1 * lengths_x[x] * lengths_z[z] * cell.D_surface[g][3]
							
						self.A[cell_id * ng + g, cell_id * ng + g] += lengths_x[x] * lengths_z[z] * cell.D_surface[g][3]
						
						#Transport to below cell
						if z != 0:
							self.A[cell_id * ng + g, (cell_id - cx * cy) * ng + g] = -1 * lengths_x[x] * lengths_y[y] * cell.D_surface[g][4]
							
						self.A[cell_id * ng + g, cell_id * ng + g] += lengths_x[x] * lengths_y[y] * cell.D_surface[g][4]
						
						#Transport to above cell
						if z != cz - 1:
							self.A[cell_id * ng + g, (cell_id + cx * cy) * ng + g] = -1 * lengths_x[x] * lengths_y[y] * cell.D_surface[g][5]
							
						self.A[cell_id * ng + g, cell_id * ng + g] += lengths_x[x] * lengths_y[y] * cell.D_surface[g][5]
							

		self.A = self.A.tocsr()
		self.a_end = time.time()

		
						
	def Compute_M(self):
 		print "Creating gain matrix..."
 		self.m_start = time.time()
 		ng = self.num_groups
 		cx = self.mesh.num_cells_x
 		cy = self.mesh.num_cells_y
 		cz = self.mesh.num_cells_z
 		
 		self.M = self.M.tolil()
 		self.M = self.M * 0.0

		for z in range(cz):
		
			for y in range(cy):
		
				for x in range(cx):
			
					cell = self.mesh.cells[z,y,x]
					cell_id = (z * cx * cy + y * cx + x)
			
					for g in range(ng):
				
						for e in range(ng):
					
							self.M[cell_id * ng + g, cell_id * ng + e] = cell.volume * cell.Material.NuSig_f[e] * cell.Material.Chi[g]

		self.M = self.M.tocsr()
		self.m_end = time.time()

#####################################
# Iterated Guesses/Solve for Fluxes # 
#####################################
	def Solve(self):
		print "Solving matrix problem..."
		self.s_start = time.time()
		ng = self.num_groups
		cx = self.mesh.num_cells_x
		cy = self.mesh.num_cells_y
		cz = self.mesh.num_cells_z
		nc = cx * cy * cz
	
		self.b_old = self.M.dot(self.Phi)
		self.b_old_sum = self.b_old.sum()
		
		
		self.b_old = self.b_old/(self.b_old_sum/(ng*nc))
		self.b_old_sum = ng * nc
		self.err = 1
		self.tol = 1e-8
		self.iteration = 0
		
		while self.err > self.tol:
 			self.Phi = spsolve(self.A,self.b_old) #Implement own written solver package
 			self.b_new = self.M.dot(self.Phi)
 			self.b_new_sum = self.b_new.sum()
 			self.b_new = self.b_new/(self.b_new_sum/(ng*nc))
 			self.keff = self.b_new_sum/self.b_old_sum
 			self.b_old = self.b_new
 			self.err = fabs(self.keff-self.keff_old)
 			self.keff_old = self.keff
 			self.iteration += 1
 			#print "Iteration: " + str(self.iteration) + "\tEffective K Value: " + str(self.keff)
 		
 		self.s_end = time.time()
 		
 		r_time = self.mesh.time
 		d_time = self.d_end - self.d_start
 		a_time = self.a_end - self.a_start
 		m_time = self.m_end - self.m_start
 		s_time = self.s_end - self.s_start
 		self.total_time = r_time + d_time + a_time + m_time + s_time
 		
 		print "Effective K Value: " + str(self.keff)
		print "Took " + str(self.iteration) + " iterations"
  		print "Took " + str(r_time) + " seconds to refine the mesh"
  		print "Took " + str(d_time) + " seconds to compute the D values for surfaces"
  		print "Took " + str(a_time) + " seconds to create the loss matrix"
  		print "Took " + str(m_time) + " seconds to create the gain matrix"
  		print "Took " + str(s_time) + " seconds to solve the matrix problem"
 		print "Took " + str(self.total_time) + " seconds to solve the entire problem"
 		
 	
 	