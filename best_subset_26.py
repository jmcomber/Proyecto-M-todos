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

###


mu = 0
for i in range(p):
	for j in range(p):
		if i != j:
			aux = sum(A[k][i] * A[k][j] for k in range(n))
			mu = max(aux, mu)

# mu[k] -> mu * k
# gamma_k -> 1 - mu * (k-1)

# M_l = (1 / (1 - mu * (k - 1))) * sum(sum(A[i][j] * b[j] for i in range(n)) for j in range(k))

# M_U = min((1 / (1 - mu * (k-1))) * (sum(sum(A[i][j] ** 2 * b[i] ** 2 for i in range(n)) for j in range(p))) ** .5, (1 / (1 - mu * (k-1))) * sum(b[i] * b[i] for i in range(n)) ** .5)


###

#### Estimación parámetros con warm start

warm_start = Model("Warm start")
warm_start.setParam("TimeLimit", 13)

beta = warm_start.addVars(range(p), vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, name="beta")

#z es al revés (1 - z) que en el modelo, para poder hacer el SOS-1
z = warm_start.addVars(range(p), vtype=GRB.BINARY, name="z")

for i in range(p):
	warm_start.addSOS(GRB.SOS_TYPE1, [beta[i], z[i]])

warm_start.addConstr(quicksum(1 - z[i] for i in range(p)) <= k)

obj = quicksum(.5 * (b[i] - quicksum(A[i][j] * beta[j] for j in range(p))) * (b[i] - quicksum(A[i][j] * beta[j] for j in range(p)))  for i in range(n))

warm_start.setObjective(obj, GRB.MINIMIZE)

warm_start.optimize()

TAU = 2

NORMA_INF = max(beta[i].X for i in range(p))
NORMA_UNO = sum(abs(beta[i].X) for i in range(p))

# M_U = min(TAU * NORMA_INF, M_U)
# M_l = min(k * M_U, M_l)
M_U = TAU * NORMA_INF
M_l = k * M_U
M_gamma_U = max(sum(sorted(A[i])[-k:]) for i in range(n)) * M_l
M_gamma_l = max(sum(abs(max(A[i][j] for j in range(p))) * M_U for i in range(n)), (k ** .5) * (sum(b[i] * b[i] for i in range(n))) ** .5)

print(M_U, M_l, M_gamma_U, M_gamma_l)


# print(M_gamma_l)
# # M_gamma_Us = [238.35401813972464, 238.85152027373076, 230.18119340244277, 248.17724886296833, 246.09344072238267, 246.65775451423895, 232.57244789152384, 245.85814447640254, 237.97132669614294, 235.26236711896883, 237.85691164556255, 240.4921284603472, 249.04395457770556, 235.71474722627948, 243.79720469174754, 245.19055216254134, 243.25298473422902, 237.957830458291, 240.76353373487447, 233.64517174896542, 238.44096008669692, 242.8566020426555, 251.30370151734064, 244.88564027583703, 245.8601366962522, 231.41321500117186, 236.73186956967504, 233.16188867305465, 238.2782980482847, 240.41971944172042, 229.5253699375951, 239.45681867750505, 250.4352559693745, 249.4139627149234, 243.99637485359466, 242.3407428461869, 245.79937090102953, 238.1983851753523, 239.7309868355601, 235.1061594463158]
k = 1
while True:
	print("CASO K = ", k)
	model = Model("Best subset - {}".format(k))
	# model.setParam("TimeLimit", 30)

	beta = model.addVars(range(p), vtype=GRB.CONTINUOUS, lb=-M_U, ub=M_U, name="beta")

	#z es al revés (1 - z) que en el modelo, para poder hacer el SOS-1
	z = model.addVars(range(p), vtype=GRB.BINARY, name="z")

	gamma = model.addVars(range(n), vtype=GRB.CONTINUOUS, lb=-M_gamma_U, ub=M_gamma_U, name="gamma")
	# gamma = model.addVars(range(n), vtype=GRB.CONTINUOUS, name="gamma")

	abs_beta = model.addVars(range(p), vtype=GRB.CONTINUOUS, lb=0, ub=M_U, name="abs_beta")

	abs_gamma = model.addVars(range(n), vtype=GRB.CONTINUOUS, lb=0, ub=M_gamma_U, name="abs_gamma")
	# abs_gamma = model.addVars(range(n), vtype=GRB.CONTINUOUS, lb=0, name="abs_gamma")


	for i in range(p):
		model.addSOS(GRB.SOS_TYPE1, [beta[i], z[i]])

	model.addConstr(quicksum(1 - z[i] for i in range(p)) <= k)

	model.addConstr(sum(abs_beta[i] for i in range(p)) <= M_l)

	model.addConstrs((gamma[i] == quicksum(A[i][j] * beta[j] for j in range(p)) for i in range(n)), name="")

	model.addConstr(quicksum(abs_gamma[i] for i in range(n)) <= M_gamma_l)

	model.addConstrs((abs_beta[i] >= beta[i] for i in range(p)), name="")
	model.addConstrs((abs_beta[i] >= -beta[i] for i in range(p)), name="")

	model.addConstrs((abs_gamma[i] >= gamma[i] for i in range(n)), name="")
	model.addConstrs((abs_gamma[i] >= -gamma[i] for i in range(n)), name="")

	# model.addConstrs((gamma[i] <= M_gamma_Us[i] for i in range(n)), name="")

	# acc_A = [sum(A[i][j] for i in range(n)) for j in range(p)]
	# root_k = k ** .5
	# norm_b = sum(b[i] * b[i] for i in range(n)) **.5

	# model.addConstrs((quicksum(abs_gamma[i] for i in range(n)) - acc_A[j] \
	# 	* quicksum(abs_beta[r] for r in range(p)) <= 0 for j in range(p)), name="")

	# model.addConstr(quicksum(abs_gamma[i] for i in range(n)) <= root_k * norm_b)

	obj = .5 * quicksum(gamma[i] * gamma[i] for i in range(n)) - \
		quicksum(quicksum(A[i][j] * b[i] for i in range(n)) * beta[j] for j in range(p)) + \
		.5 * quicksum(b[i] * b[i] for i in range(n))

	model.setObjective(obj, GRB.MINIMIZE)

	model.optimize()
	if model.objval <= 0.001:
		break
	k += 1

# for v in model.getVars():
# 	if v.X != 0 and "beta" in v.varName and "abs" not in v.varName:
# 		print("{} {}".format(v.varName, v.X))


