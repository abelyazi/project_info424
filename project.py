from __future__ import division
import re
from stat import FILE_ATTRIBUTE_ARCHIVE
import pyomo.environ as pyo
from pyomo.environ import *



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

"""Récuperation des données à input dans le modèle"""
FILE = './Instances/bin_pack_20_2.dat'
data.load(filename=FILE)

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

"""Affichage des résultats"""

instance = model.create_instance(data)
opt = pyo.SolverFactory('glpk')
result = opt.solve(instance)

print(result)

"""Affichage des variable x pour lesqueslles on a un résultat"""

for i in instance.x:
    if pyo.value(instance.x[i]) > 0:
        print(instance.x[i], "de valeur : ", pyo.value(instance.x[i]))