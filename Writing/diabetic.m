%% Generalised Smooth of Insulin data
clear;
clc;

%% Load the Data
x = (0:.5:4.5)';
%normal
%y  = [70,150,165,145,90,75,65,75,80,75]';
%diabetic
y  = [100,185,210,220,195,175,105,100,85,90]';

n=length(y);

scatter(x,y)

%% Set up the basis with a knot at every point

knots    = x;
norder   = 4;
nbasis   = n + (norder - 2);

rng      = [x(1),x(end)];
basisobj = create_bspline_basis(rng, nbasis, norder, knots);
B        = eval_basis(x, basisobj);

%% Indicate the derivatives to be estimated
% [D0,D1]
Lcoef        = [1,1]';

%% Generalised Smoothoing
argvals=x;  Term=Lcoef;
tic;
[ beta_opt, C_opt, y_hat, var_beta, pars] = Gen_Pen_Non(argvals, y, basisobj, Term);
toc

%Optimal beta
display('Estimated Parameters of the ODE:');
display(beta_opt);

%Variance of beta
display('Estimated Variance of the Parameters of the ODE:');
display((var_beta));

display('Estimated Standard Deviation of the Parameters of the ODE:');
display(sqrt(var_beta));

%Optimal lambda
display('Estimated complexity Parameter, Sum of Squared Errors and Degrees of Freedom');
display(pars);

%Plot
figure();
plot(x,y,'*k')
hold on;
plot(x,y_hat,'-k')