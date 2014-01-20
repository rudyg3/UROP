# All are of similar geometries to top-down view, but the bottom-most layer (0-20)
# is all water (4), whereas the top-most layer (360-380) alternates between
# water with control rod (5) and water (4), where the control rod is located between
# x = 30 and x = 50 and y = 30 and y = 50. For the layer (280-360), everything is the same
# as the top-down view, but there is a control rod (3) again between x = 30 and x = 50 
# and y = 30 and y = 50. 

import os
import sys
sys.path.insert(0, os.getcwd() + '/../src')
import numpy as np
import matplotlib.pyplot as plt
import Nuclear as nuke
import Mesh
from Solver import Solver