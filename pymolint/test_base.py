import numpy as np 
from itertools import permutations
from pymolint.base import *


def test_base():
    """
    """
    find_contacts([(1,1,1),(2,2,2),(3,3,3)],[1,2,3],[(1.5,1.5,1.5),(2.5,2.5,2.5),(3.5,3.5,3.5)],[4,5,6])
