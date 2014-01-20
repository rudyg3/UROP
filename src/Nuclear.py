import numpy as np
from pprint import pprint
from math import fabs
from math import sqrt
from time import sleep


class Materials(object):
	"""Contains data regarding the materials, separated into arrays for different groups"""
	def __init__(self, num_groups, ID):
		self.ID = ID
		self.Chi = np.zeros(num_groups)
		self.D = np.zeros(num_groups)
		self.Sig_a = np.zeros(num_groups)
		self.NuSig_f = np.zeros(num_groups)
		self.Sig_s = np.zeros((num_groups,num_groups))
		self.num_groups = num_groups
		
		
	def Set_Chi(self, Chi):
		try:
			assert isinstance(Chi, np.ndarray)
		except:
			print "Chi must be array"
		self.Chi = Chi
		
	def Set_D(self, D):
		try:
			assert isinstance(D, np.ndarray)
		except:
			print "D must be array"
		self.D = D
		
	def Set_Sig_a(self, Sig_a):
		try:
			assert isinstance(Sig_a, np.ndarray)
		except:
			print "Sig_a must be array"
		self.Sig_a = Sig_a
		
	def Set_NuSig_f(self, NuSig_f):
		try:
			assert isinstance(NuSig_f, np.ndarray)
		except:
			print "NuSig_f must be array"
		self.NuSig_f = NuSig_f
	
	def Set_Sig_s(self, Sig_s):
		try:
			assert isinstance(Sig_s, np.ndarray)
		except:
			print "Sig_s must be array"
		self.Sig_s = Sig_s
		
	def Set_ID(self, ID):
		self.ID = ID
	
class Cell(object):
	"""Contains data regarding the cell"""
	
	def __init__(self, num_groups):
		self.D_surface = np.zeros((num_groups, 6))
		self.Flux = np.zeros((num_groups,1))
		self.num_groups = num_groups
		self.volume = 0
		
	def Set_Material(self, Material):
		"""Material: a pointer to an instance of the materials class"""
		try:
			assert isinstance(Material, Materials)
		except:
			"Material must be instance of Materials class"
		self.Material = Material
	