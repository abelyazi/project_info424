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
                    self.get_dad().set_ub(ub_min)
                    new_bound = True

                if lb_min != dad.get_lb():
                    self.get_dad().set_lb(lb_min)
                    new_bound = True

                if new_bound==True:
                    self.get_dad().update_lb_ub()


"""Etablissement fonction objective"""
def obj_expression(m):
    return pyo.summation(m.y)

"""Etablissement des contraintes"""
def x_constraint_rule(m, p):
    return sum(m.x[p,b] for b in m.B) == 1

def sx_constraint_rule(m, b):
    return sum(m.size[p] * m.x[p,b] for p in m.P) <= m.cap*m.y[b]

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
def solve_bp_lp(instance_name):
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
    model.OBJ = pyo.Objective(rule=obj_expression)
    model.xpbConstraint = pyo.Constraint(model.P, rule=x_constraint_rule)
    model.sxConstraint = pyo.Constraint(model.B, rule=sx_constraint_rule)
    file = './Instances/' + instance_name
    data.load(filename=file)

    instance = model.create_instance(data)
    opt = pyo.SolverFactory('glpk')

    start_time = time.time()
    opt.solve(instance)
    tts = time.time() - start_time
    print("LP relaxation of " + instance_name + " solved in %s seconds." % tts)

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

    return obj, x, y



def solve_instnc_for_BnB(instnc):
    try:
        instance = instnc
        opt = pyo.SolverFactory('glpk')
        #start_time = time.time()
        opt.solve(instance)
        #print("LP solved in %s seconds." % (time.time() - start_time))
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


def branch_and_bound(instance_name, branching_scheme, valid_inequalities,time_limit):
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
    model.OBJ = pyo.Objective(rule=obj_expression)
    model.xpbConstraint = pyo.Constraint(model.P, rule=x_constraint_rule)
    model.sxConstraint = pyo.Constraint(model.B, rule=sx_constraint_rule)

    file = './Instances/' + instance_name
    data.load(filename=file)
    instance = model.create_instance(data)
    
    for j in range(len(list(instance.size))-1):
        instance.Constraints.add(instance.y[j] >= instance.y[j+1] )
    

    if valid_inequalities !=0:
        instance_cp = deepcopy(instance)
        sol_cp = solve_instnc_for_BnB(instance_cp)
        x_cp = sol_cp[1]
        x_cp_new = []
        for b in range(len(x_cp)):
            valeur = 0
            for p in range(len(x_cp)):
                if (x_cp[p][b] > 0) and (x_cp[p][b] < 1):
                    if x_cp[p][b] >= valeur:
                        valeur = x_cp[p][b]
                        u = p, b
            x_cp_new.append(u)

        for elem_cp in x_cp_new:
            index_p = elem_cp[0]
            index_b = elem_cp[1]
            expr = 0
            s_cp = pyo.value(instance.size[index_p])
            for k in range(len(x_cp)):
                new_coef = (pyo.value(instance.size[k])) / s_cp
                expr += floor(new_coef) * instance.x[(k, index_b)]
            coef_y = floor((pyo.value(instance.cap)) / s_cp)
            instance.Constraints.add(expr <= coef_y * instance.y[index_b])
    
    q = deque([])
    root_node = Node(None,None,None,0,None,None,None)
    q.append(root_node)
    start_time1 = time.time()
    while (len(q) != 0) and (time.time() < start_time1 + time_limit):
        if branching_scheme==0:
            current = q.pop()  # Depth
        elif branching_scheme==1:
            current = q.popleft() # Breadth

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
        current_sol = solve_instnc_for_BnB(current_instance)
        lb = ceil(deepcopy(current_sol[0])) #résultat de la fonction objective arrondis à la valeur supérieur car on traitre des solutions entières
        current.set_lb(lb)
        #print("THIS IS LB ",lb)

        
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
            #print("THIS IS UB",ub)
        
        # step 4 if solution trouvée possède soit ub!=lb, soit sol faisable non entière
        
        ## step 4.1 sélection de la variable sur laquelle on va imposer la constrainte >= 1 et <= 0
        ## pour les noeuds fils
            
            current.update_lb_ub()
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
                                    
                    fils_constraints=[[(a,b),0],[(a,b),1]]
                    ## step 4.2 création des deux noeuds et on les ajoute à la queue q
                    current.set_right_child(Node(None,None,fils_constraints[0],current.get_level()+1,current,None,None))
                    current.set_left_child(Node(None,None,fils_constraints[1],current.get_level()+1,current,None,None))
                    q.append(current.get_right_child())
                    q.append(current.get_left_child())
                    if current.get_level()!=0:
                        #print(current.get_dad().get_lb())
                        #print(current.get_dad().get_ub())
                        current.update_lb_ub()
                        #print(current.get_dad().get_lb())
                        #print(current.get_dad().get_ub())

                        temp_dad = deepcopy(current)
                        bool_val = False
                        while temp_dad != None:
                            #print("this is the lb and ub of node of lvl : ",temp_dad.get_level()," ",temp_dad.get_lb()," ",temp_dad.get_ub())
                            if temp_dad.get_lb() == temp_dad.get_ub():
                                bool_val = True
                            temp_dad = temp_dad.get_dad()
                        if bool_val == True:
                            buffer = q.pop()
                            buffer = q.pop()
    
    while current.get_dad() != None:
        current = current.get_dad()
    return current.get_lb(), current.get_ub()



instance_name="bin_pack_100_1.dat"
start_time0 = time.time()
l_b0, u_b0 = branch_and_bound(instance_name,0,0,60)
time_BnB0 = time.time() - start_time0

start_time1 = time.time()
l_b1, u_b1 = branch_and_bound(instance_name,1,0,60)
time_BnB1 = time.time() - start_time1

print("Depth first: ")
print(time_BnB0, ' seconds')
print("Solution lb et ub du root node",l_b0, u_b0)

print('\n')

print("Breadth first: ")
print(time_BnB1, ' seconds')
print("Solution lb et ub du root node",l_b1, u_b1)