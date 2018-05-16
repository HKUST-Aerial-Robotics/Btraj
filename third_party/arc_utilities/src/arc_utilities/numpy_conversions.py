import numpy as np
from geometry_msgs.msg import *

def ListPointsToNpArray(points):
    arr = np.empty(shape=(3, len(points)))
    for ind in range(len(points)):
        arr[:,ind] = [points[ind].x, points[ind].y, points[ind].z]
    return arr