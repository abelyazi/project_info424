from __future__ import division
from collections import deque
from pickle import TRUE

import re
from stat import FILE_ATTRIBUTE_ARCHIVE
import time
from typing import MutableMapping, MutableSequence, MutableSet
import pyomo.environ as pyo
from pyomo.environ import *

class Node:
    def __init__(self,ub,lb,cnstr,lvl,dad):
      self.uper_b = ub
      self.lower_b = lb
      self.constraint=cnstr
      self.tree_lvl= lvl 
      self.parent_node = dad

    def get_level(self):
        return self.level
    
    def get_dad(self):
        return self.parent_node
    
    def get_constraint(self):
        return self.constraint

    def get_ub(self):
        return self.uper_b

    def get_lb(self):
        return self.lower_b

    def set_ub(self,ub):
        self.uper_b = ub

    def set_lb(self,lb):
        self.lower_b,lb
    


"""Etablissement variables et paramètres"""
model = pyo.AbstractModel()
data = pyo.DataPortal(model=model)
model.I = pyo.Set()
model.P = pyo.Set(initialize=model.I)
model.B = pyo.Set(initialize=model.I)
model.size = pyo.Param(model.P)
model.cap = pyo.Param()
model.x = pyo.Var(model.P,model.B, domain=pyo.NonNegativeReals, bounds=(0,1))
model.y = pyo.Var(model.B, domain=pyo.NonNegativeReals, bounds=(0,1))
model.Constraints = pyo.ConstraintList()


"""Etablissement fonction objective"""
def obj_expression(m):
    return pyo.summation(m.y)
model.OBJ = pyo.Objective(rule=obj_expression)

"""Etablissement des contraintes"""
def x_constraint_rule(m, p):
    return sum(m.x[p,b] for b in m.B) == 1
model.xpbConstraint = pyo.Constraint(model.B, rule=x_constraint_rule)

def sx_constraint_rule(m, b):
    return sum(m.size[p] * m.x[p,b] for p in m.P) <= m.cap*m.y[b]
model.sxConstraint = pyo.Constraint(model.P, rule=sx_constraint_rule)

def check_results_x(instnc):
    """Affichage des variable x pour lesqueslles on a un résultat > 0"""
    for i in instnc.x:
        if (pyo.value(instnc.x[i]) > 0.001):
            print(instnc.x[i], "de valeur : ", pyo.value(instnc.x[i]))
def check_results_y(instnc):
    """Affichage des boites utilisées"""
    for i in instnc.y:
        if (pyo.value(instnc.y[i]) > 0):
            print(instnc.y[i], "=", pyo.value(instnc.y[i]))

def check_x_is_real(instnc):
    for i in instnc.x:
        if (pyo.value(instance.x[i]) > 0.001) and (pyo.value(instance.x[i]) < 1) :
            return True
    return False

def check_y_is_real(instnc):
    for i in instnc.y:
        if (pyo.value(instance.y[i]) > 0.001) and (pyo.value(instance.y[i]) < 1) :
            return True
    return False    

# Résolution d'une instance quelconque
def solve_instnc(instnc):
    instance = instnc
    opt = pyo.SolverFactory('glpk')
    start_time = time.time()
    opt.solve(instance)
    print("LP solved in %s seconds." % (time.time() - start_time))
    x = []
    for p in range(len(list(instance.size))):
        x.append([])
        for b in range(len(list(instance.size))):
            x[p].append(0)
    for p in range(len(list(instance.size))):
        for b in range(len(list(instance.size))):
            x[p][b] = pyo.value(instance.x[p, b])
    y = []
    for i in instance.y:
        y.append(pyo.value(instance.y[i]))

    obj = instance.OBJ()
    #solution = obj,x,y
    return obj, x, y

# Résolution
def solve_bp_lp(instance_name):
    instance_name="bin_pack_20_2.dat"
    file = './Instances/' + instance_name
    data.load(filename=file)

    instance = model.create_instance(data)
    opt = pyo.SolverFactory('glpk')
    start_time = time.time()
    opt.solve(instance)
    print("LP solved in %s seconds." % (time.time() - start_time))


    x = []
    for p in range(len(list(instance.size))):
        x.append([])
        for b in range(len(list(instance.size))):
            x[p].append(0)
    for p in range(len(list(instance.size))):
        for b in range(len(list(instance.size))):
            x[p][b] = pyo.value(instance.x[p, b])

    y = []
    for i in instance.y:
        y.append(pyo.value(instance.y[i]))

    obj = instance.OBJ()
    #solution = obj,x,y
    return obj, x, y, instance



# Résolution du noeud racine
solution = solve_bp_lp("bin_pack_20_2.dat")
instance = solution[3]
lvl = 0
visited = []
q = deque([])
q.append(Node(None,None,None,lvl,None))
while (len(q) != 0)and(lvl<450):
    current = q.popleft()  # Breadth
    current_instance = instance

    # step 1 on remplit la liste de containtes avec les contraintes des noeuds parents jusqu'au noeud racine
    constraints=[]
    temp_current = current
    while temp_current.get_constraint() != None:
        constraints.append(temp_current.get_constraint())
        if temp_current.get_dad() != None:
            temp_current = temp_current.get_dad()
    #print(constraints)

    # step 2 établissement lb : résolution du problem en ajoutant la liste constraintes au problème de base
    for i in constraints:
        expr = 0
        expr += current_instance.x[i[0]]
        if i[1] == 0:
            current_instance.x[i[0]].fix(0)
        elif i[1] == 1:
            current_instance.x[i[0]].fix(1)
    
    current_sol = solve_instnc(current_instance)
    lb = current_sol[0] #résultat de la fonction objective
    current.set_lb(lb)
    #print(lb)

    # step 3 etablissement ub : réparation de la solution trouvé et établissement du lb

    # step 4 if solution trouvée possède soit ub!=lb, soit sol faisable non entière
    
    ## step 4.1 sélection de la variable sur laquelle on va imposer la constrainte >= 1 et <= 0
    ## pour les noeuds fils
    x = current_sol[1]
    temp = 1
    for i in range(len(x)):
        for j in range(len(x[i])):
            val = x[i][j]
            if (val>0.01) and (val<1):
                val2 = min(val,1-val)
                if val2<temp:
                    temp = val2
                    a,b = i,j
                    #print(current_instance.x[(a,b)])
    print(x[a][b]) 
    fils_constraints=[[(a,b),1],[(a,b),0]]

    ## step 4.2 création des deux noeuds et on les ajoute à la queue q
    node1= Node(None,None,fils_constraints[0],lvl+1,current)
    node2= Node(None,None,fils_constraints[1],lvl+1,current)
    q.append(node1)
    q.append(node2)

    # step 5 
    visited.append(current)
    
    lvl+=1
    

   

"""
if current[0] == 7 or current[0] == 6:
    print(f"Found Goal {current[0]} with cost: {current[1]}")
    break
"""   



"""Réparation des variables x
while check_x_is_real(instance):
    temp=1
    nb = 0
    key = 0
    for i in instance.x:
        val = pyo.value(instance.x[i])
        if (val>0.01) and (val<1):
            val2 = min(val,1-val)
            if val2<temp:
                temp = val2
                nb = val
                key = i
                print(key)
                print(val)
    instance.x[key].fix(round(nb))
    results = opt.solve(instance)
    check_results_x(instance)
"""
    
"""Réparation des variables y
while check_y_is_real(instance):
    temp=1
    nb = 0
    key = 0
    for i in instance.y:
        val = pyo.value(instance.y[i])
        if (val>0) and (val<1):
            val2 = min(val,1-val)
            if val2<temp:
                temp = val2
                nb = val
                key = i
                print(key)
                print(val)
    instance.y[key].fix(round(nb))
    results = opt.solve(instance)
    check_results_x(instance)
    check_results_y(instance)
"""


"""Méthode heuristique : First Fit Decreasing
def FirstFit(I, cap, size):
    bins = [cap]
    for i in range(len(size)):
        iFit = False
        for j in range(len(bins)):
            if size[i] <= bins[j]:
                bins[j] = bins[j] - size[i]
                break
            if j == len(bins)-1:
                iFit = True

        if iFit is True:
            bins.append(cap)
            bins[len(bins)-1] -= size[i]

    return len(bins)

def FF_decreasing_heuristic(I, cap, size):
    sizeSorted = size.copy()
    sizeSorted.sort(reverse = True)
    return FirstFit(I, cap, sizeSorted)

n = 60
c = 150
w = []
for i in instance.size:
    w.append(instance.size[i])

print(FF_decreasing_heuristic(n,c,w)) 
"""