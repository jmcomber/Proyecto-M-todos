function [optval,sol,solarray] = AFOM(A,b,tau,NITER)

cpuini = cputime;

[m,n] = size(A)

xk = zeros(n,1);
xk1 = zeros(n,1);
xk2 = zeros(n,1);
zig = 0;
thetak = 1;
zk = xk;
angle = 0;

R = sqrt(norm((A*A')^-1))*norm(b)

L = tau*sqrt(n) + norm(A*A')*R + norm(A'*b)

optval = tau*norm(xk,1) + 0.5*norm(A*xk-b)**2;
B = A'*A;
bb = A'*b;

for t = 1:NITER,
  yk = (1-thetak)*xk + thetak*zk;
  g = B*yk - bb;
  for j = 1:n,
     x1 = zk(j) -(1/(L*thetak))*(g(j)+tau);
     x2 = zk(j) -(1/(L*thetak))*(g(j)-tau);
     x3 = 0;
     val1 = g(j)*x1 + tau*abs(x1)+(L*thetak/2)*(x1-zk(j))*(x1-zk(j));
     val2 = g(j)*x2 + tau*abs(x2)+(L*thetak/2)*(x2-zk(j))*(x2-zk(j));
     val3 = g(j)*x3 + tau*abs(x3)+(L*thetak/2)*(x3-zk(j))*(x3-zk(j));
     if val1 <= min(val2,val3) 
        zk(j) = x1;
     else 
        if val2 <= min(val1,val3)
           zk(j) = x2;
        else 
           zk(j) = x3;
        endif
     endif
  end
  zk;
  xk2 = xk1;
  xk1 = xk;
  xk = (1-thetak)*xk + thetak*zk;
  thetak = 2/(1+sqrt(1+4/(thetak*thetak)));
  aux = norm(A*xk-b);
  optval = tau*norm(xk,1) + 0.5*aux**2;
  error_residual_relativo = max(abs(A*xk-b))/max(abs(b));
 # if t >= 3
 #   zig = norm(xk2-xk1)/norm(xk-xk2);
 #   angle = ((xk2-xk1)'*(xk-xk1)/(norm(xk2-xk1)*norm(xk-xk1)));
 # endif
  solarray(t,1) = t;
  solarray(t,2) = optval;
  solarray(t,3) = norm(xk,1);
  solarray(t,4) = error_residual_relativo;
  solarray(t,5) = zig;
  solarray(t,6) = angle;
  # printf("%6i  %10.6f %10.6f %10.6f %10.6f \n",t,optval,thetak,norm(xk,1),error_residual_relativo);
  fflush(stdout);
  end
 

cpufinal = cputime-cpuini
sol = xk;


endfunction