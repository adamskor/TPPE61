%% Heston

nu0   = 0.1  ; 
kappa = 10   ; 
eta_   = 0.1  ; 
theta_ = 1.0  ; 
rho   = -0.80; 
%************************************************
paramS = [nu0;kappa;eta;theta;rho];

S0  = 100;
K   = 100;
TTM = 1;
r   = 0;
div = 0;
modellnamn = 'Heston';



paramS_ = [nu0;kappa;eta_;theta_;rho]


[pOrg, grad_] = heston_([100,110]', [100, 115]', [0, 0]', [0, 0]', [1.2, 1.2]', paramS_)
% [p, grad_] = heston_TestGrad(S0, K, r, div, TTM, paramS_)
