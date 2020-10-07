clc
clear
close all

d = 50;

Theta1 = randn(d,d)+ 50*eye(d);

Theta1 = max(Theta1,Theta1');
Sigma1 = inv(Theta1);

Mask = [rand(d,d)>0.95];

Mask = double(Mask);

Delta = randn(d,d);

Delta = Delta.*Mask;

Delta = max(Delta,Delta');

for i = 1:d
    Delta(i,i) = 0;
end

Theta2 = Theta1 + Delta;
Sigma2 = inv(Theta2);

subplot(2,3,1)
imagesc(Theta1)
hold on
subplot(2,3,2)
imagesc(Theta2)
subplot(2,3,3)
imagesc(Delta)

%generate data
n = 200;
mu = zeros(d,1);
X1 = mvnrnd(mu,Sigma1,n);
hSigma1 = cov(X1);
norm(hSigma1-Sigma1)

X2 = mvnrnd(mu,Sigma2,n);

hSigma2 = cov(X2);
norm(hSigma2-Sigma2)

lambda = 0.004;
[hDelta] = differential_graph(hSigma1,hSigma2,lambda);

norm(hDelta-Delta)

subplot(2,3,6)
imagesc(hDelta)
