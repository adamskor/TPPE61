%%%Enligt Lourakis 2005   
S0 = readmatrix('optionData.xlsx','Range','B2:B2');
K = readmatrix('optionData.xlsx','Range','I2:I37');
TTM = readmatrix('optionData.xlsx','Range','K2:K37');
C_star = readmatrix('optionData.xlsx','Range','H2:H37');

%[f, p] = gridSearch(C_star, S0, K , TTM)

nSamples = 20;
nu0LB = 0.05; nu0UB = 0.95; nu0 = (nu0UB-nu0LB).*rand(nSamples,1) + nu0LB;
kappaLB = 0.5; kappaUB = 5; kappa = (kappaUB-kappaLB).*rand(nSamples,1) + kappaLB;
etaLB = 0.05; etaUB = 0.95; eta = (etaUB-etaLB).*rand(nSamples,1) + etaLB;
thetaLB = 0.05; thetaUB = 0.95; theta = (thetaUB-thetaLB).*rand(nSamples,1) + thetaLB;
rhoLB = -0.1; rhoUB = -0.9; rho = (rhoUB-rhoLB).*rand(nSamples,1) + rhoLB;
figure(1)
check = zeros(5,nSamples);
objValues = zeros(1,nSamples);
e_pVec = zeros(numel(K),nSamples);

for i = 1:nSamples
    p0 = [nu0(i); kappa(i); eta(i); theta(i); rho(i)];
    [p, f, fIter, e_p] = LevenbergFuncGeneral(@func2, p0, C_star, S0, K , TTM);
    check(:,i) = p;
    objValues(i) = f;
    e_pVec(:,i) = e_p;
end

scatter(1:1:nSamples ,objValues)


