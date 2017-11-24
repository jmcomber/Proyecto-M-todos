load Grupo6.mat

NITER = 2000;
lambdas = [0.1, 0.5, 0.7];

for p = 1:3
  tau = lambdas(p);
  printf("\n\nCaso lambda = %6i\n", tau)
  [optval,sol,solarray] = FISTA(A,b,tau,NITER);
  umbrales = [0.7, 0.8, 0.9, 1.0, 1.1];
  for j = 1:5
    cont_vars = 300;
    for i = 1:300
      if sol(i) <= umbrales(j)
        sol(i) = 0.0;
        cont_vars -= 1;
      endif
    endfor
    nuevo_optimo = 0.0;
    for h = 1:39
      aux = 0;
      otro = 0;
      for k =1:300
        aux += A(h*300 + k) * sol(k);
      endfor
      aux -= b(h);
      aux = aux ^ 2;
      otro += aux;
      aux += tau * abs(b(h));
      nuevo_optimo += aux;
    endfor
    printf("Desajuste es %6i\n", otro / nuevo_optimo)
    printf("Con umbral %6i, %6i variables distintas de cero: valor optimo: %6i\n", umbrales(j), cont_vars, nuevo_optimo)  
  endfor
endfor