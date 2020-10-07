function [Theta] = differential_graph(Sigma1,Sigma2,lambda)

Q = kron(Sigma1,Sigma2);

d = size(Sigma1,1);

b = Sigma1 - Sigma2;

b = b(:);

nVars = d*d;

w_init = zeros(nVars,1);

funObj = @(w)DGLoss(w,Q,b);

params = [];
params.verbose = 0;
[w] = L1GeneralCompositeGradientAccelerated(funObj,w_init,lambda,params);

Theta = reshape(w,[d,d]);

Theta_T = Theta';
ind = abs(Theta) < abs(Theta_T);
Theta(ind) = Theta_T(ind);

end

function [f,g] = DGLoss(w,Q,b)

f = 1/2 *w'*Q*w -  w'*b;

g = Q*w - b;

end
