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
    def __init__(self,ub,lb,cnstr,lvl,dad,l_child,r_child):
      self.uper_b = ub
      self.lower_b = lb
      self.constraint=cnstr
      self.tree_lvl= lvl 
      self.parent_node = dad
      self.left_child = l_child
      self.right_child = r_child
      

    def get_level(self):
        return self.tree_lvl
    
    def get_dad(self):
        return self.parent_node
    
    def get_constraint(self):
        return self.constraint

    def get_ub(self):
        return self.uper_b
    
    def get_left_child(self):
        return self.left_child
    
    def get_right_child(self):
        return self.right_child
    
    def get_lb(self):
        return self.lower_b

    def set_dad(self,dad):
        self.parent_node = dad
    
    def set_left_child(self,l_child):
        self.left_child=l_child
    
    def set_right_child(self,r_child):
        self.right_child=r_child

    def set_ub(self,ub):
        self.uper_b = ub

    def set_lb(self,lb):
        self.lower_b = lb
    
    def update_lb_ub(self):
        if self.get_dad() != None:
            dad = self.get_dad()
            left_child = dad.get_left_child()
            right_child = dad.get_right_child()
            if (right_child.get_ub() != None) and (left_child.get_ub() != None):
                ub_min = min(right_child.get_ub(),left_child.get_ub())
                lb_min = min(right_child.get_lb(),left_child.get_lb())

                if (right_child.get_lb()==0) or (left_child.get_lb()==0):
                    lb_min = right_child.get_lb() + left_child.get_lb()
                
                new_bound = False

                if ub_min != dad.get_ub():
                    print("chgt ub")
                    self.get_dad().set_ub(ub_min)
                    new_bound = True

                if lb_min != dad.get_lb():
                    print("chgt lb")
                    self.get_dad().set_lb(lb_min)
                    new_bound = True

                if new_bound==True:
                    self.get_dad().update_lb_ub()

            
    


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

def check_results_x(instnc):
    """Affichage des variable x pour lesqueslles on a un résultat > 0"""
    for i in instnc.x:
        if (pyo.value(instnc.x[i]) > 0):
            print(instnc.x[i], "de valeur : ", pyo.value(instnc.x[i]))
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

def check_dad_ub_are_bigger_than_lb_current(node):
    current_n = deepcopy(node)
    l_b = node.get_lb()
    while current_n != None :
        #print("Lower Bound du noeud au lvl", node.get_level(),":", l_b)
        #print("Upper Bound du noeud au lvl", current_n.get_level(),":", current_n.get_ub())
        if l_b > current_n.get_ub():
            return False
        current_n = deepcopy(current_n.get_dad())
    return True
            


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

instance_name="bin_pack_90_0.dat"
file = './Instances/' + instance_name
data.load(filename=file)
instance = model.create_instance(data)

#ajout des contraintes yj >= yj+1
for j in range(len(list(instance.size))-1):
    instance.Constraints.add(instance.y[j] >= instance.y[j+1] )


# Résolution du noeud racine
iteration = 0
visited = []
q = deque([])
tree_node_list = deque([])
root_node = Node(None,None,None,0,None,None,None)
tree_node_list.append(root_node)
q.append(root_node)
start_time1 = time.time()
while (len(q) != 0) and (iteration<20000) and (time.time() < start_time1 + 600):
    current = q.pop()  # Depth
    #current = q.popleft() # Breadth
    
    """
    if current == root_node:
        tree_node = root_node
    # code pour mettre 
    elif current != root_node:
        if tree_node.get_level() < current.get_level():
            temp_node = current.get_dad()
            while (temp_node.get_left_child()!=current) and (temp_node.get_right_child()!=current):
                temp_node = current.get_dad()
            if temp_node.get_left_child()==current:
                print("trouvey1")
                tree_node = tree_node.get_left_child()
            if temp_node.get_right_child()==current:
                print("trouvey2")
                tree_node = tree_node.get_right_child()

        elif tree_node.get_level() >= current.get_level():
            tree_node = tree_node.get_dad()
            while (tree_node.get_left_child()!=current) and (tree_node.get_right_child()!=current):
                tree_node = tree_node.get_dad()
            if tree_node.get_left_child()==current:
                print("trouvey3")
                tree_node = tree_node.get_left_child()
            if tree_node.get_right_child()==current:
                print("trouvey4")
                tree_node = tree_node.get_right_child()
    """    


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
        if const[1] == 0:
                current_instance.Constraints.add( current_instance.x[key] <= 0 )
        elif const[1] == 1:
                current_instance.Constraints.add( current_instance.x[key] >= 1 )


    # step 2 établissement lb : résolution du problem en ajoutant la liste constraintes au problème de base
    current_sol = solve_instnc(current_instance)
    lb = ceil(deepcopy(current_sol[0])) #résultat de la fonction objective arrondis à la valeur supérieur car on traitre des solutions entières
    
    current.set_lb(lb)

    #tree_node.set_lb(lb)

    print("THIS IS LB ",lb)

    
    # step 3 etablissement ub : réparation de la solution trouvé et établissement du lb
    if (lb!=0):
        x = deepcopy(current_sol[1])
        y = deepcopy(current_sol[2])
        ## step 3.1 on multiplie chaque élement xij par la taille du paquet i
        x_temp = deepcopy(x)
        x_ij = []
        for p in range(len(x_temp)):
            for b in range(len(x_temp[p])):
                if x_temp[p][b] > 0: 
                    x_temp[p][b] *= pyo.value(instance.size[p])
                    x_ij.append((x_temp[p][b],(p,b)))

            
        
        # step 3.2 on parcourt les éléments de x_ij (triés) et on ajoute les paquets dans les boites une par une
        x_ij = sorted(x_ij,reverse=True)
        #print(x_ij)
        capa_y = [0 for i in range(len(y))]
        presence_x = [0 for i in range(len(x))]
        capa = pyo.value(current_instance.cap)
        boxes_used = 0
        for elem in x_ij:
            packet_index = elem[1][0]
            box_index = 0
            if presence_x[packet_index]==0:
                presence_x[packet_index]=1
                while ((current_instance.size[packet_index] + capa_y[box_index]) > capa) and box_index < len(y):
                    box_index += 1
                if capa_y[box_index]==0:
                    boxes_used +=1
                capa_y[box_index] += current_instance.size[packet_index]

        ub = boxes_used
        current.set_ub(ub)
        #tree_node.set_ub(ub)
        print("THIS IS UB",ub)
    
    # step 4 if solution trouvée possède soit ub!=lb, soit sol faisable non entière
    
    ## step 4.1 sélection de la variable sur laquelle on va imposer la constrainte >= 1 et <= 0
    ## pour les noeuds fils

        if ub>lb and (check_dad_ub_are_bigger_than_lb_current(current)==True):
            if check_x_is_real(current_instance)==True:            
                temp=0
                for i in range(len(x)):
                    for j in range(len(x[i])):
                        val = x[i][j]
                        if (val>0) and (val<1):
                            if val>temp:
                                temp = val
                                a,b = i,j
                                #print(current_instance.x[(a,b)])
                print((a,b))               
                print(x[a][b]) 
                fils_constraints=[[(a,b),0],[(a,b),1]]
                ## step 4.2 création des deux noeuds et on les ajoute à la queue q
                #node_right = Node(None,None,fils_constraints[0],current.get_level()+1,current,None,None)
                #node_left = Node(None,None,fils_constraints[1],current.get_level()+1,current,None,None)
                current.set_right_child(Node(None,None,fils_constraints[0],current.get_level()+1,current,None,None))
                current.set_left_child(Node(None,None,fils_constraints[1],current.get_level()+1,current,None,None))
                #tree_node.set_right_child(node_right)
                #tree_node.set_left_child(node_left)
                q.append(current.get_right_child())
                q.append(current.get_left_child())
                if current.get_level()!=0:
                    print(current.get_dad().get_lb())
                    print(current.get_dad().get_ub())
                    current.update_lb_ub()
                    print(current.get_dad().get_lb())
                    print(current.get_dad().get_ub())

                    temp_dad = deepcopy(current)
                    bool_val = False
                    while temp_dad != None:
                        print("this is the lb and ub of node of lvl : ",temp_dad.get_level()," ",temp_dad.get_lb()," ",temp_dad.get_ub())
                        if temp_dad.get_lb() == temp_dad.get_ub():
                            bool_val = True
                        temp_dad = temp_dad.get_dad() 
                    """"
                    while  (temp_dad.get_dad() != None):
                        print(1)
                        if temp_dad.get_ub() == temp_dad.get_lb():
                            print(2)
                            bool_val = True
                            temp_dad = temp_dad.get_dad()
                    """
                    if bool_val == True:
                        current = q.pop()
                        current = q.pop()

                
                """
                if tree_node.get_level()!=0:
                    print(tree_node.get_dad().get_lb())
                    print(tree_node.get_dad().get_ub())
                    tree_node.update_lb_ub()
                    print(tree_node.get_dad().get_lb())
                    print(tree_node.get_dad().get_ub()) 
                """
                #current = update_lb_ub(current)
                #check_results_x(current_instance)
                
    iteration+=1    


    # step 5 

 
    #visited.append(current)
    
print("SOLUTION")
while current != None:
    print("this is the lb and ub of node of lvl : ",current.get_level()," ",current.get_lb()," ",current.get_ub())
    current = current.get_dad() 
#for node in visited:
#    print(node.get_level())
#check_results_x(current_instance)
#check_results_y(current_instance)
#print(visited[0].get_ub())
#print(visited[0].get_lb())
   

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