from __future__ import division
import pyomo.environ as pyo

model = pyo.AbstractModel()



model.I = pyo.Set()

model.P = pyo.Set(initialize=model.I)
model.B = pyo.Set(initialize=model.I)

model.size = pyo.Param(model.P)
model.cap = pyo.Param(domain=pyo.NonNegativeIntegers)

# the next line declares variable x indexed by the set P and B
model.x = pyo.Var(model.P,model.I, domain=pyo.NonNegativeReals, bounds=(0,1))

# the next line declares variable y indexed by the set B
model.y = pyo.Var(model.I, domain=pyo.NonNegativeReals, bounds=(0,1))

def obj_expression(m):
    return pyo.summation(m.y)

model.OBJ = pyo.Objective(rule=obj_expression)

def x_constraint_rule(m, p):
    # return the expression for the constraint for i
    return sum(m.x[p,b] for b in m.B) == 1

# the next line creates one constraint for each member of the set model.I
model.xpbConstraint = pyo.Constraint(model.B, rule=x_constraint_rule)

def sx_constraint_rule(m, b):
    # return the expression for the constraint for i
    return sum(m.size[p] * m.x[p,b] for p in m.P) <= m.cap*m.y[b]

# the next line creates one constraint for each member of the set model.I
model.sxConstraint = pyo.Constraint(model.P, rule=sx_constraint_rule)

