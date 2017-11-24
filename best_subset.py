from gurobipy import *


n, p = 40, 300
k = 5

file = open("Grupo6.mat", "r")
A, b = [], []
j = 0
for i in file.readlines():
	if j >= 5 and j < 45:
		line = [float(k) for k in i.strip().split(" ")]
		A.append(line)
	elif j > 50 and j < 91:
		b.append(float(i))
	j += 1


M_U = 137.21707657
M_l = 32798.8048


## ESTIMACIÓN PARÁMETROS (INICIO)

# u_mas = []
# u_menos = []

# # model.objval

# UB = .5 * quicksum((b[i] - sum(A[i][j] for j in range(k))) * (b[i] - sum(A[i][j] for j in range(k)))  for i in range(n))

# for i in range(p):
# 	print("ESTIMACIÓN {}".format(i))
# 	m_mas = Model("u_mas_iter")
# 	m_menos = Model("u_menos_iter")

# 	beta_mas = m_mas.addVars(range(p), vtype=GRB.CONTINUOUS, name="beta_mas")
# 	beta_menos = m_menos.addVars(range(p), vtype=GRB.CONTINUOUS, name="beta_menos")
	
# 	m_menos.addConstr(quicksum((b[i] - sum(A[i][j] * beta_menos[j] for j in range(p))) * (b[i] - sum(A[i][j] * beta_menos[j] for j in range(p)))  for i in range(n)) <= 2 * UB)
# 	m_mas.addConstr(quicksum((b[i] - sum(A[i][j] * beta_mas[j] for j in range(p))) * (b[i] - sum(A[i][j] * beta_mas[j] for j in range(p)))  for i in range(n)) <= 2 * UB)
	
# 	m_mas.setObjective(beta_mas[i], GRB.MAXIMIZE)
# 	m_menos.setObjective(beta_menos[i], GRB.MINIMIZE)
	
# 	m_mas.optimize()
# 	u_mas.append(m_mas.objval)
# 	m_menos.optimize()
# 	u_menos.append(m_menos.objval)

# M_Us = [max(abs(u_mas[i]), abs(u_menos[i])) for i in range(p)]
# M_U = max(M_Us)
# M_l = sum(M_Us)

# print("M_U = {}".format(M_U))
# print("M_l = {}".format(M_l))


## ESTIMACIÓN PARÁMETROS (FIN)

model = Model("Best subset |k|")

beta = model.addVars(range(p), vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, name="beta")

#z es al revés (1 - z) que en el modelo, para poder hacer el SOS-1
z = model.addVars(range(p), vtype=GRB.BINARY, name="z")

abs_beta = model.addVars(range(p), vtype=GRB.CONTINUOUS, lb=0, ub=M_U, name="abs_beta")

for i in range(p):
	model.addSOS(GRB.SOS_TYPE1, [beta[i], z[i]])

model.addConstr(quicksum(1 - z[i] for i in range(p)) <= k)
# model.addConstr(quicksum(z[i] for i in range(p)) >= p - k)

model.addConstrs((-M_U <= beta[i] for i in range(p)), name="extra (1)")
model.addConstrs((M_U >= beta[i] for i in range(p)), name="extra (2)")

model.addConstr(quicksum(abs_beta[i] for i in range(p)) <= M_l)

model.addConstrs((abs_beta[i] >= beta[i] for i in range(p)), name="")
model.addConstrs((abs_beta[i] >= -beta[i] for i in range(p)), name="")

# obj = quicksum(.5 * (b[i] - quicksum(A[i][j] * beta[j] for j in range(p))) * (b[i] - quicksum(A[i][j] * beta[j] for j in range(p)))  for i in range(n))

obj = .5 * quicksum(quicksum(A[h][j] * beta[j] for j in range(p)) * quicksum(A[h][j] * beta[j] for j in range(p)) for h in range(n)) - \
	quicksum(quicksum(A[i][j] * b[i] for i in range(n)) * beta[j] for j in range(p)) + \
	.5 * quicksum(b[i] * b[i] for i in range(n))


model.setObjective(obj, GRB.MINIMIZE)

model.optimize()

# for v in model.getVars():
# 	if v.X != 0:
# 		print("{} {}".format(v.varName, v.X))


