function [theta, f] = LevenbergFunc(options, theta_0)
    theta = theta_0;

        
modellnamn = 'Heston';

m = 5;
n = size(options, 1);



eps1 = 10e-15;
eps2 = 10e-15;
eps3 = 10e-15;

%%[p_,grad] = mexOption_ps2('Heston', S0, K, r, div, TTM, theta_0);
J = calcJacobian(options, theta, m, n);
r = calcR(options, theta, n);
my = 10^-3*max(diag(J));
nu = 2;
I = eye(m);

[theta, f] = Leven(J, r, theta, my, nu, I, options, n, m, eps1, eps2, eps3);
function [theta, f] = Leven(J, r, theta, my, nu, I, options, n ,m, eps1, eps2, eps3)
    while true
        [deltaF, deltaL, theta, thetaStep] = compute(theta, r, J, my, I, n, options);
        if deltaL > 0 && deltaF > 0
            J = calcJacobian(options, theta, m, n);
            r = calcR(options, theta, n);
        else
            my = my*nu;
            nu = 2*nu;
        end
        con1 = calcNormR(r);
        con2 = max(sum(abs(J'*r)));
        con3 = norm(thetaStep)/norm(theta);
        if (con1 <= eps1)
            break;
        end
        if (con2 <= eps2)
            break;
        end
        if (con3 <= eps3)
            break;
        end
        

    end
    f = calcF(r);
end

function [deltaF, deltaL, theta, thetaStep] = compute(theta, r, J, my, I, n, options)
    thetaStep = (J*J' + my*I)\calcGradF(J, r);
    theta = theta + thetaStep;
    rNext = calcR(options, theta, n);
    deltaL = thetaStep'*(my*thetaStep + J*r);
    deltaF = calcNormR(r) - calcNormR(rNext);
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
end