from __future__ import division
from cmath import inf
from collections import deque
from copy import deepcopy
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
        return self.tree_lvl
    
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
model.xpbConstraint = pyo.Constraint(model.P, rule=x_constraint_rule)

def sx_constraint_rule(m, b):
    return sum(m.size[p] * m.x[p,b] for p in m.P) <= m.cap*m.y[b]
model.sxConstraint = pyo.Constraint(model.B, rule=sx_constraint_rule)


def check_results_x(matrix_x):
    """Affichage des variable x pour lesqueslles on a un résultat > 0"""
    for i in range(len(matrix_x)):
        for j in range(len(matrix_x[i])):
            if (matrix_x[i][j] > 0):
                print("x[",i,"][",j,"] de valeur : ", matrix_x[i][j])

def check_results_y(instnc):
    """Affichage des boites utilisées"""
    for i in instnc.y:
        if (pyo.value(instnc.y[i]) > 0):
            print(instnc.y[i], "=", pyo.value(instnc.y[i]))

def check_x_is_real(instnc):
    for i in instnc.x:
        if (pyo.value(instnc.x[i]) > 0) and (pyo.value(instnc.x[i]) < 1) :
            return True
    return False

def check_y_is_real(instnc):
    for i in instnc.y:
        if (pyo.value(instnc.y[i]) > 0) and (pyo.value(instnc.y[i]) < 1) :
            return True
    return False    

# Résolution d'une instance quelconque
def solve_instnc(instnc):
    try:
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
    except:
        print("Solution Infeasible")
        obj = (0,0,0)
        return obj 



instance_name="bin_pack_50_0.dat"
file = './Instances/' + instance_name
data.load(filename=file)
instance = model.create_instance(data)

#ajout des contraintes yj >= yj+1
for j in range(len(list(instance.size))-1):
    instance.Constraints.add(instance.y[j] >= instance.y[j+1] )

solution0 = solve_instnc(instance)
print(solution0[0])
x = solution0[1]
y = solution0[2]

x_temp = deepcopy(x)

for p in range(len(x_temp)):
    for b in range(len(x_temp[p])):
        x_temp[p][b] = x_temp[p][b] * pyo.value(instance.size[p])


print('\n')


x_ij = []
for p in range(len(x_temp)):
    for b in range(len(x_temp[p])):
        if x_temp[p][b] != 0:
            x_ij.append((x_temp[p][b],(p,b)))

#print(sorted(x_ij,reverse=True))


capa_y = [0 for i in range(len(y))]
presence_x = [0 for i in range(len(x))]
capa = 150
for elem in x_ij:
    packet_index = elem[1][0]
    box_index = 0
    if presence_x[packet_index]==0:
        presence_x[packet_index]=1
        while ((instance.size[packet_index] + capa_y[box_index]) > capa) and box_index < len(y):
            box_index += 1
        capa_y[box_index] += instance.size[packet_index]

print(capa_y)
print(presence_x)


     
    




"""
# Résolution du noeud racine
iteration = 0
visited = []
q = deque([])
q.append(Node(None,None,None,0,None))
while (len(q) != 0)and(iteration<1000):
    current = q.pop()  # Depth
    #current = q.popleft() # Breadth
    current_instance = deepcopy(instance)


    # step 1 on ajoute la contrainte du current_node à l'instance du noeud

    constraints=[]
    temp_current = deepcopy(current)

    while temp_current.get_constraint() != None:
        constraints.append(deepcopy(temp_current.get_constraint()))
        if temp_current.get_dad() != None:
            temp_current = deepcopy(temp_current.get_dad())

    for const in constraints:
        key = const[0] #(a,b) 
        if key[0] == 20:
            if const[1] == 0:
                    current_instance.Constraints.add( current_instance.y[key[1]] <= 0 )
            elif const[1] == 1:
                    current_instance.Constraints.add( current_instance.y[key[1]] >= 1 )
        elif key[0] != 20:
            if const[1] == 0:
                    current_instance.Constraints.add( current_instance.x[key] <= 0 )
            elif const[1] == 1:
                    current_instance.Constraints.add( current_instance.x[key] >= 1 )

       
        

    # step 2 établissement lb : résolution du problem en ajoutant la liste constraintes au problème de base
    current_sol = solve_instnc(current_instance)
    lb = deepcopy(current_sol[0]) #résultat de la fonction objective
    current.set_lb(lb)
    print(lb)

    # step 3 etablissement ub : réparation de la solution trouvé et établissement du lb

    # step 4 if solution trouvée possède soit ub!=lb, soit sol faisable non entière
    
    ## step 4.1 sélection de la variable sur laquelle on va imposer la constrainte >= 1 et <= 0
    ## pour les noeuds fils
    xy_matrix=[]
    if lb != 0:
        if check_x_is_real(current_instance)==True or (check_y_is_real(current_instance)==True):
            xy_matrix = deepcopy(current_sol[1])
            y = deepcopy(current_sol[2])
            xy_matrix.append(y)
            temp=0
            for i in range(len(xy_matrix)):
                for j in range(len(xy_matrix[i])):
                    val = xy_matrix[i][j]
                    if (val>0) and (val<1):
                        if val>temp:
                            temp = val
                            a,b = i,j
                            #print(current_instance.x[(a,b)])
            print((a,b))               
            print(xy_matrix[a][b]) 
            fils_constraints=[[(a,b),0],[(a,b),1]]
            ## step 4.2 création des deux noeuds et on les ajoute à la queue q
            q.append(Node(None,None,fils_constraints[0],current.get_level()+1,deepcopy(current)))
            q.append(Node(None,None,fils_constraints[1],current.get_level()+1,deepcopy(current)))
        


    # step 5 
    visited.append(deepcopy(current))
    
    iteration+=1
#for node in visited:
#    print(node.get_level())

#check_results_x(current_instance)
#check_results_y(current_instance)
"""
   

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