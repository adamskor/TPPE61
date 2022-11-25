


%           S0  K   TTM  r  div marketPrice
 options = [100 90  0.2  0  0   3;
            100 95  0.4  0  0   5;
            100 100 0.6  0  0   12;
            100 105 1.5  0  0   16;
            100 110 0.6  0  0   23;
            ];
        
modellnamn = 'Heston';

m = 5;
n = size(options, 1);

theta = [0.1 10 0.1 1.0 -0.80]';

eps1 = 10e-9;
eps2 = 10e-9;
eps3 = 10e-9;

%%[p_,grad] = mexOption_ps2('Heston', S0, K, r, div, TTM, theta_0);
J = calcJacobian(options, theta, m, n);
r = calcR(options, theta, n);
my = max(diag(J));
nu = 2;
I = eye(m);

[theta] = Leven(J, r, theta, my, nu, I, options, n, m, eps1, eps2, eps3);
function [theta] = Leven(J, r, theta, my, nu, I, options, n ,m, eps1, eps2, eps3)
    %J = J; r = r; my = my; nu = nu; options = options; n = n; m = m;
    %I = I; eps1 = eps1; eps2 = eps2; eps3 = eps3;
    while true
        thetaStep = (J*J' + my*I)\calcGradF(J, r);
        theta = theta + thetaStep;
        rNext = calcR(options, theta, n);
        normR = calcNormR(r);
        deltaL = thetaStep'*(my*thetaStep + J*r);
        deltaF = calcNormR(r) - calcNormR(rNext);
        if deltaL > 0 && deltaF > 0
            J = calcJacobian(options, theta, m, n);
            r = rNext;
        else
            my = my*nu;
            nu = 2*nu;
            r = rNext;
            Leven(J, r, theta, my, nu, I, options, n ,m, eps1, eps2, eps3);
        end
        con1 = normR;
        con2 = max(sum(abs(J')));
        con3 = norm(thetaStep)/norm(theta);
        if (con1 < eps1) || (con2 < eps2) || (con3 < eps3)
            break;
        end
        

    end

end




function J = calcJacobian(options, theta, m, n) 
    J = zeros(m, n);
    for i = 1:n
        [~ ,grad] = mexOption_ps2('Heston', options(i, 1), options(i, 2), options(i, 4), options(i, 5), options(i, 3), theta);
        J(:, i) = grad;
    end
end

function r = calcR(options, theta, n)
    r = zeros(n, 1);
    for i = 1:n
        [C ,~] = mexOption_ps2('Heston', options(i, 1), options(i, 2), options(i, 4), options(i, 5), options(i, 3), theta);
        r(i) = C - options(i, 6);
    end
end

function normR = calcNormR(r)
    normR = norm(r);
end


function gradF = calcGradF(J, r)
    gradF = J*r;
end

function f = calcF(r)
    f = 0.5*r'*r;
end
