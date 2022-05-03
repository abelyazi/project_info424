from collections import deque


class Node:
   def __init__(self, instnc):
      self.uper_b = None
      self.lower_b = None
      self.new_constraint=None
      self.tree_lvl= None

def branchAndBound(instance_name, branching_scheme,valid_inequalities,time_limit):

    
    q = deque([])
    while len(q) != 0:
        current = q.popleft()  # Breadth
        # current = q.pop()  # Depth
        if current[0] == 7 or current[0] == 6:
            print(f"Found Goal {current[0]} with cost: {current[1]}")
            break
        for neighbors in graph[current[0]]:
            q.append((neighbors[0], current[1] + neighbors[1]))