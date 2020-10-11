'''
Package for aranging miscellaneous things
'''
__author__ = "Giyeok Lee"
__email__ = "lgy4230@yonsei.ac.kr"
__date__ = "Sep 14, 2020"
__maintainer__ = "Giyeok Lee"
__version__ = "2020.09.14"
__copyright__ = "Giyeok Lee"

from numpy import cos, sin, tan, sqrt
import numpy as np
import glob
import os
import argparse
import copy

# class AR: @staticmethod
def length(X):
    from math import sqrt
    '''
    calculate the length of vector. Regardless of the direction.
    X : np.array, dtype=float
    Output (Results) : float type
    '''
    X=np.array(X,dtype='d')
    temp = 0
    for i in X:
        temp+=i**2
    return sqrt(temp)

def species(X, overall = True):
    '''
    For comparing just species and its' orders (not number of atoms), we use not np.array, but list.
    [input]
    X : list, tuple, np.array any iterable type.
    overall : don't care about sequences. True if you jlust want to know species
    [output]
    species : list type.
    '''
    species=[]
    for i in X:
        if len(species)==0:
            species.append(i)
        else:
            if overall:
                if i not in species:
                    species.append(i)
            if not overall:
                if i!=species[-1]:
                    species.append(i)
    return species

def npsort(X, axis=-1):
    '''
    Maintaining the coordinates of ions, and sort along axis.
    X = np.array type (n x 3 matrix)
    '''
    temp=[]
    for i in np.argsort(X[:,axis]):
        temp.append(X[i])
    return np.array(temp)
        
def tolsort(X, tol = 1e-5):
    '''
    When X have values with small difference,
    and if you want to remove that values, using this function

    for now, only 1D array is considered
    '''
    before_tol = np.sort(X)
    after_tol = before_tol[[True]+list(np.abs(np.diff(before_tol))>tol)]
    return after_tol        
    
def hor_translate(X, shift = 0):
    A = X[shift:]
    B = X[:shift]
    return np.hstack((X[shift:], X[:shift]))    

def dir2car(coord, cell_vec):
    '''Direct to Cartesian
    [input] : coord, cell_vec
    |-> coord : coordination. @np.array(dtype="d")
    |-> cell_vec : cell vector. @np.array(dtype='d')
    [output] : new coordination @np.array(dtype='d')
    '''
    return np.matmul(coord, cell_vec)

def car2dir(coord, cell_vec):
    '''Cartesian to Direct
    [input] : coord, cell_vec
    |-> coord : coordination. @np.array(dtype='d')
    |-> cell_vec : cell vector. @np.array(dtype='d')
    [output] : new coordination @np.array(dtype='d')
    '''
    inv_cell_vec = np.linalg.inv(cell_vec)
    return np.matmul(coord, inv_cell_vec)
  
def to_new_cell(model, decimals = 13):
    new_cell = np.zeros((3,3))
    a, b, c, alpha, beta, gamma = model.get_cell_lengths_and_angles()
    alpha, beta, gamma = np.array([alpha, beta, gamma]) * np.pi / 180
    new_cell[0,0] = a
    new_cell[1] = b * cos(gamma), b * sin(gamma), 0
    new_cell[2,0] = c * cos(beta)
    new_cell[2,1] = c * (cos(alpha) - cos(beta) * cos(gamma)) / (sin(gamma))
    new_cell[2,2] = sqrt(c ** 2 - new_cell[2,0] ** 2 - new_cell[2,1] ** 2)
    if decimals != None:
        new_cell = new_cell.round(decimals = decimals)
    tmp = model.copy()
    tmp.set_cell(new_cell)
    tmp.set_scaled_positions(model.get_scaled_positions())
    return tmp

def pop(model, index = [-1]):
    '''
    [Description]
    : Remove selected atoms from model(ase.Atoms instance)
      and return popped atoms with original unit cell.

    [args]
    * index : list of integer

    [Output]
    * return the structure after popping
    '''
    tmp = model.copy()
    removed = []
    for i in model:
        if i.tag in index:
            removed.append(i.index)
    del tmp[removed]
    return tmp

def dirlist():
    all = glob.glob("*")
    direc = [x for x in all if os.path.isdir(x)]
    return direc