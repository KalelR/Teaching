%% Phytoplankton model with limiting resources
% 
%
%%
clear ;  close all ; clc
%% Parameters
% n = k = 3
S = [10, 10, 10];
K = [1.0 0.75 0.25 0.7 0.2 0.65 0.68 0.38 0.46;
     0.25 1.0 0.75 0.2 1.01 0.55 0.83 1.10 0.85;
     0.75 0.25 1 1.10 0.7 0.95 0.6 0.5 0.77];

C = [0.10 0.20 0.15 0.05 0.01 0.40 0.30 0.20 0.25;
    0.15 0.10 0.20 0.15 0.30 0.35 0.25 0.02 0.35;
     0.20 0.15 0.10 0.25 0.05 0.20 0.40 0.15 0.10];
K = K(1:3, 1:3);
C = C(1:3, 1:3);


% n = k = 5
S = [6, 10, 14, 4, 9];
K = [0.39 0.34 0.30 0.24 0.23 0.41 0.20 0.45 0.14 0.15 0.38 0.28;
     0.22 0.39 0.34 0.30 0.27 0.16 0.15 0.05 0.38 0.29 0.37 0.31;
     0.27 0.22 0.39 0.34 0.30 0.07 0.11 0.05 0.38 0.41 0.24 0.25;
     0.30 0.24 0.22 0.39 0.34 0.28 0.12 0.13 0.27 0.33 0.04 0.41;
     0.34 0.30 0.22 0.20 0.39 0.40 0.50 0.26 0.12 0.29 0.09 0.16];

C = [0.04 0.04 0.07 0.04 0.04 0.22 0.10 0.08 0.02 0.17 0.25 0.03;
     0.08 0.08 0.08 0.10 0.08 0.14 0.22 0.04 0.18 0.06 0.20 0.04;
     0.10 0.10 0.10 0.10 0.14 0.22 0.24 0.12 0.03 0.24 0.17 0.01;
     0.05 0.03 0.03 0.03 0.03 0.09 0.07 0.06 0.03 0.03 0.11 0.05;
     0.07 0.09 0.07 0.07 0.07 0.05 0.24 0.05 0.08 0.10 0.02 0.04];
K = K(:, 1:5);
C = C(:, 1:5);


[k, n] = size(C);
r = zeros(n, 1); r(:) = 1;
D = 0.25 ;
m = zeros(n, 1); m(:) = D;
%p = [K, C, S, r, m, D];
%# p = Dict(:K=>K, :C=>C, :S=>S, :r=>r, :m=>m, :D=>D)

% ics 
R0 = S;
N0 = 0.1 + (1:n) ./ 100;
u0 = [N0 R0];



%% Simulation
time_sim    = linspace(0, 300, 1000); 
[t,states] = ode45(@(t,states) resourcecompetition(t,states, K, C, S, r, m, D),time_sim, u0);

hold on
for i=1:n 
    plot(t, states(:, i)) 
end
%plot3(states(:,1), states(:,2), states(:,3))

function du = resourcecompetition(~,states, K, C, S, r, m, D)
    [k, n] = size(C);
    N = states(1:n);
    R = states(n+1:n+k);
    for i=1:n
        du(i,1) = N(i)*(mu(i, R,r,K) - m(i) );
    end
    for j=1:k
        du(n+j, 1) = D*(S(j) - R(j)) - resource_consumption(j, C, R, r, K, N);
    end
end

function sum = resource_consumption(j, C, R, r, K, N)
    sum = 0.0;
    n = length(N);
    for i=1:n
        sum = sum + C(j,i) * mu(i, R, r, K) * N(i);
    end
end

function mu = mu(i, R, r, K) 
    mu = min(arrayfun(@monod_equation, r, R, K(:,i)));
end
function mu = monod_equation(r, R, Ks)
    mu = r*R/(Ks+R);
end