import numpy as np                                                                                                               import scipy.spatial as ss                                                                                                       from itertools import combinations
from itertools import permutations

from base import *


def test_base():
    """
    """
    find_contacts([(1,1,1),(2,2,2),(3,3,3)],[1,2,3],[(1.5,1.5,1.5),(2.5,2.5,2.5),(3.5,3.5,3.5)],[4,5,6])