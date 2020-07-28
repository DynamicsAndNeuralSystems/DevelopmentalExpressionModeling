% Try to optimize rho, d0scale


%-------------------------------------------------------------------------------
% SET INITIAL CONDITIONS
%-------------------------------------------------------------------------------

% -------Optimize for just lambda and A:
% d0Scaling_0 = 7.51; % d0Scaling_0 = 3.84;
% rho_0 = 0.13; % rho_0 = 0.21;

% -------Optimize for lambda, A, and f0:
d0Scaling_0 = 8.40; % d0Scaling_0 = 3.84;
rho_0 = 0.16; % rho_0 = 0.21;

X0 = [d0Scaling_0,rho_0];

fprintf(1,'Down-hill search\n');
% Down-hill search from a given starting point:
options = struct('Display','iter','MaxIter',50,'TolFun',1e-1,'TolX',1e-1,'PlotFcns',@optimplotfval);
errorFun = @(x)computeError(x(1),x(2));
paramOpt = fminsearch(errorFun,X0,options);

% We get: 15.5798,0.4164
