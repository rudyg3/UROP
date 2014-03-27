import numpy as np
import matplotlib.pyplot as plt
import Nuclear as nuke
import time
from math import *


class Mesh(object):
	def __init__(self, lengths_x, lengths_y, lengths_z, num_groups):
		self.num_cells_x = len(lengths_x)
		self.num_cells_y = len(lengths_y)
		self.num_cells_z = len(lengths_z)
		self.lengths_x = lengths_x
		self.lengths_y = lengths_y
		self.lengths_z = lengths_z
		self.boundary = ["REFLECTIVE"] * 6
		self.time = 0
		self.num_cells = self.num_cells_x * self.num_cells_y * self.num_cells_z
		self.num_groups = num_groups
		self.cells = [nuke.Cell(num_groups) for i in range(self.num_cells_x * self.num_cells_y * self.num_cells_z)]
		self.cells = np.reshape(self.cells,(self.num_cells_z, self.num_cells_y, self.num_cells_x))
		
		for z in range(self.num_cells_z):
			for y in range(self.num_cells_y):
				for x in range(self.num_cells_x):
					self.cells[z][y][x].volume = lengths_x[x] * lengths_y[y] * lengths_z[z]
		
	def refineMesh(self, cell_size, dimensions=['x','y','z']):
		print "\n\n\nRefining mesh..."
		r_start = time.time()
		if 'x' in dimensions:
			lengths_x = []
			for i in self.lengths_x:
				new_length_x = int(floor(i/cell_size))
				new_lengths_x = []
				
				for j in range(new_length_x):
					new_lengths_x.append(cell_size)
				
				for k in new_lengths_x:
					lengths_x.append(k)
					
		else:
			lengths_x = self.lengths_x
					
		if 'y' in dimensions:
			lengths_y = []
			for i in self.lengths_y:
				new_length_y = int(floor(i/cell_size))
				new_lengths_y = []
				
				for j in range(new_length_y):
					new_lengths_y.append(cell_size)
				
				for k in new_lengths_y:
					lengths_y.append(k)	
					
		else:
			lengths_y = self.lengths_y
		
		if 'z' in dimensions:
			lengths_z = []
			for i in self.lengths_z:
				new_length_z = int(floor(i/cell_size))
				new_lengths_z = []
				
				for j in range(new_length_z):
					new_lengths_z.append(cell_size)
					
				for k in new_lengths_z:
					lengths_z.append(k)
				
		else:
			lengths_z = self.lengths_z
					

		mesh = Mesh(lengths_x, lengths_y, lengths_z, self.num_groups)
		
		for z in range(mesh.num_cells_z):
		
			for y in range(mesh.num_cells_y):
			
				for x in range(mesh.num_cells_x):
			
					xval = (sum(mesh.lengths_x[0:x]) + sum(mesh.lengths_x[0:x+1]))/2
					yval = (sum(mesh.lengths_y[0:y]) + sum(mesh.lengths_y[0:y+1]))/2
					zval = (sum(mesh.lengths_z[0:z]) + sum(mesh.lengths_z[0:z+1]))/2
				
					for xn in range(self.num_cells_x):
						if (xval >= sum(self.lengths_x[0:xn]) and xval <= sum(self.lengths_x[0:xn+1])):
							break
				
					for yn in range(self.num_cells_y):
						if (yval >= sum(self.lengths_y[0:yn]) and yval <= sum(self.lengths_y[0:yn+1])):
							break
						
					for zn in range(self.num_cells_z):
						if (zval >= sum(self.lengths_z[0:zn]) and zval <= sum(self.lengths_z[0:zn+1])):
							break
			
					mesh.cells[z][y][x].Set_Material(self.cells[zn][yn][xn].Material)
		
		r_end = time.time()
		mesh.time = r_end - r_start
		return mesh
			
	def Get_Material_ID(self,cell_instance):
		return cell_instance.Material.ID
	
	def Set_Boundaries(self, boundaries):
		self.boundary = boundaries