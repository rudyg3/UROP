from TwoLRA import LRA
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np

Go = LRA()

mesh_refinement = np.ones(14)
times = np.ones(14)

#for j in range(1,11):
mesh_refinement_prime = []
times_prime = []
for i in range(1,15):
	mesh_refinement_prime.append(float(i))
	times_prime.append(Go.Compute(float(i)))
	# for k in range(len(mesh_refinement_prime)):
# 		mesh_refinement[k] += mesh_refinement_prime[k]
# 		times[k] += times_prime[k]
# 		
# mesh_refinement = mesh_refinement/10.0
# times = times/10.0


plt.close()
#plt.xkcd()
plt.xlabel("Mesh refinement (Centimeters)")
plt.ylabel("Time taken to solve LRA problem (Seconds)")
plt.title("Efficiency of 2D LRA Solver")
plt.semilogy(mesh_refinement_prime,times_prime, 'bo', linestyle='none')
ax = plt.gca()
ax.set_xticks([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0])
ax.yaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
plt.savefig("/Users/rudygarcia/Desktop/School/UROP/3D-diffusion/LRA/Effiency.pdf")
plt.show()

#Fix the disparity in iteration 1, Make sure method of taking averages works (this means optimizing the matrix solve method)