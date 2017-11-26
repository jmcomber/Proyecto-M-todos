import numpy as np


def grad(beta):
	# retorna vector del gradiente de norma A * beta - b al cuadrado evaluado en beta
	vector = []
	ajuste = np.array(A*beta - b)
	for i in range(p):
		col_i = np.squeeze(np.asarray(A[:, i]))
		ajuste = np.squeeze(np.asarray(ajuste))
		vector.append(2 * np.dot(col_i, ajuste))
	return vector

n, p = 40, 300


N_ITER = 10

file = open("Grupo6.mat", "r")
A, b = [], []
j = 0
for i in file.readlines():
	if j >= 5 and j < 45:
		line = [float(m) for m in i.strip().split(" ")]
		A.append(line)
	elif j > 50 and j < 91:
		b.append(float(i))
	j += 1


A = np.matrix(A)
b = np.transpose(np.matrix(b))

R = (np.linalg.norm(np.linalg.inv(A * np.transpose(A)))**.5) * np.linalg.norm(b, ord=2)
L = np.linalg.norm(A * np.transpose(A)) * R + np.linalg.norm(np.transpose(A)*b, ord=2)

print("R", R)
print("L", L, "\n")

results = []
for k in [5, 10, 15, 20, 100, 290]:
	#beta_0
	beta = np.transpose(np.matrix([1 if i < k else 0 for i in range(p)]))

	#no hay valor previo para calcular epsilon al principio
	prev = float('inf')

	actual = np.linalg.norm(A * beta - b, ord=2) ** 2
	iteracion = 1
	while abs(actual - prev) > 10**-10: 
		gradiente = np.transpose(np.matrix(grad(beta)))
		c = beta - (1/L) * gradiente
		c = c.tolist()
		copia_ord = sorted(c, key=lambda x: -abs(x[0]))
		beta = np.matrix([c[i] if copia_ord.index(c[i]) < k else [0] for i in range(len(c))])
		prev = actual
		actual = np.linalg.norm(A * beta - b, ord=2) ** 2
		# print("Iteración {}: ajuste de {}, epsilon de {}".format(iteracion, actual, abs(actual - prev)))
		iteracion += 1
		if iteracion > N_ITER:
			results.append("k = {}: ajuste de {}, epsilon de {}".format(k, actual, abs(actual - prev)))
			break

for result in results:
	print(result)
# beta = beta.tolist()
# for i in beta:
# 	if i[0] != 0.0:
# 		print(beta.index(i), i[0])
