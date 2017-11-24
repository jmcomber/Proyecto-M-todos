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
k = 5

N_ITER = 10

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


A = np.matrix(A)
b = np.transpose(np.matrix(b))

R = (np.linalg.norm(np.linalg.inv(np.transpose(A) * A))**.5) * np.linalg.norm(b, ord=2)
L = np.linalg.norm(np.transpose(A) * A) * R + np.linalg.norm(np.transpose(A)*b, ord=2)

#beta_0
beta = np.transpose(np.matrix([1000 if i < k else 0 for i in range(p)]))

#no hay valor previo para calcular epsilon al principio
prev = float('inf')

actual = np.linalg.norm(A * beta - b, ord=2) ** 2
iteracion = 1
while abs(actual - prev) > 10**-10: 
	gradiente = np.transpose(np.matrix(grad(beta)))
	c = beta - (1/L**.5) * gradiente
	c = c.tolist()
	copia_ord = sorted(c, key=lambda x: abs(x[0]))
	nuevo_beta = [copia_ord[i][0] if copia_ord.index(copia_ord[i]) < k else 0 for i in range(len(beta))]
	copia = [0 for _ in range(len(c))]
	for i in c:
		if i[0] in nuevo_beta:
			copia[c.index(i)] = i[0]
		else:
			copia[c.index(i)] = 0
	# copia tiene los betas del nuevo beta, pero en las posiciones correctas
	beta = np.transpose(np.matrix(copia))
	prev = actual
	actual = np.linalg.norm(A * beta - b, ord=2) ** 2
	print("Iteración {}: ajuste de {}, epsilon de {}".format(iteracion, actual, abs(actual - prev)))
	iteracion += 1
	if iteracion > N_ITER:
		break

beta = beta.tolist()
for i in beta:
	if i[0] != 0.0:
		print(beta.index(i), i[0])
