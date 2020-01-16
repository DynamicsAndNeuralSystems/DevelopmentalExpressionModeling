% Try to optimize rho, d0scale

d0Scaling_0 = 3.66;
rho_0 = 0.21;
X0 = [d0Scaling_0,rho_0];

fprintf(1,'Down-hill search\n');
% Down-hill search from a given starting point:
options = struct('Display','iter','MaxIter',10,'TolFun',1e-1,'TolX',1e-1,'PlotFcns',@optimplotfval);
errorFun = @(x)computeError(x(1),x(2));
paramOpt = fminsearch(errorFun,X0,options);
