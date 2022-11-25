options = [100 90  0.2  0  0   3;
            100 95  0.4  0  0   5;
            100 100 0.6  0  0   12;
            100 105 1.5  0  0   16;
            100 110 0.6  0  0   23;
            ];
        
        
        
%[nu0;kappa;eta;theta;rho]
L = 5;
nu0_ = linspace(0.05, 0.06, L);
kappa_ = linspace(0.50, 0.51, L);
eta_ = linspace(0.05, 0.06, L);
theta_ = linspace(0.05, 0.06, L);
rho_ = linspace(-0.11, -0.1, L);
n = 1000;
nums = rand(5,n);
result = zeros(5, n);
fVec = zeros(1, n);
for i = 1:n
    [theta, f] = LevenbergFunc(options, nums(:, i));
    result(:, i) = theta;
    fVec(i) = f;
end
subplot(2, 3, 1)
plot(result(1, :))
title('nu')
subplot(2, 3, 2)
plot(result(2, :))
title('kappa')
subplot(2, 3, 3)
plot(result(3, :))
title('eta')
subplot(2, 3, 4)
plot(result(4, :))
title('theta')
subplot(2, 3, 5)
plot(result(5, :))
title('rho')
subplot(2, 3, 6)
plot(fVec)
title('f')

count = 1;
% result = zeros(5, L^5);
% for it1 = 1:L
%     for it2 = 1:L
%         for it3 = 1:L
%             for it4 = 1:L
%                 for it5 = 1:L
%                     theta = [nu0_(it1); kappa_(it2); eta_(it3); theta_(it4); rho_(it5)];
%                     result(:, count) = LevenbergFunc(options, theta);
%                     count = count + 1;
%                 end
%             end
%         end
%     end
% end

theta = [0.12 10 0.1 1.0 -0.8]';
theta2 = [0.1 10 0.1 1.0 -0.80]';
%[theta] = LevenbergFunc(options, theta)
%[theta] = LevenbergFunc(options, theta2)