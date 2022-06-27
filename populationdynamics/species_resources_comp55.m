function species_resources_comp55(r, m, D, K, C, S, T, deltat)

clf

% time-derivatives 

% mu = specific_growth_rate(r, K); % define the specific growth rates
% 
% Ndots = dNdt(mu, m); % define the species time derivatives
% 
% Rdots = dRdt(D, S, C, mu); % define the resources time derivatives
% 
% %% explicit from here on
% 
% N1dot = @(N, R) Ndots{1}(N, R);
% 
% N2dot = @(N, R) Ndots{2}(N, R);
% 
% N3dot = @(N, R) Ndots{3}(N, R);
% 
% R1dot = @(N, R) Rdots{1}(N, R);
% 
% R2dot = @(N, R) Rdots{2}(N, R);
% 
% R3dot = @(N, R) Rdots{3}(N, R);

%% manually computing the derivatives

mu1 = @(R) min([r(1)*R(1)/(K(1,1)+R(1)),r(1)*R(2)/(K(2,1)+R(2)),r(1)*R(3)/(K(3,1)+R(3)),r(1)*R(4)/(K(4,1)+R(4)),r(1)*R(5)/(K(5,1)+R(5))]);
mu2 = @(R) min([r(2)*R(1)/(K(1,2)+R(1)),r(2)*R(2)/(K(2,2)+R(2)),r(2)*R(3)/(K(3,2)+R(3)),r(2)*R(4)/(K(4,2)+R(4)),r(2)*R(5)/(K(5,2)+R(5))]);
mu3 = @(R) min([r(3)*R(1)/(K(1,3)+R(1)),r(3)*R(2)/(K(2,3)+R(2)),r(3)*R(3)/(K(3,3)+R(3)),r(3)*R(4)/(K(4,3)+R(4)),r(3)*R(5)/(K(5,3)+R(5))]);
mu4 = @(R) min([r(4)*R(1)/(K(1,4)+R(1)),r(4)*R(2)/(K(2,4)+R(2)),r(4)*R(3)/(K(3,4)+R(3)),r(4)*R(4)/(K(4,4)+R(4)),r(4)*R(5)/(K(5,4)+R(5))]);
mu5 = @(R) min([r(5)*R(1)/(K(1,5)+R(1)),r(5)*R(2)/(K(2,5)+R(2)),r(5)*R(3)/(K(3,5)+R(3)),r(5)*R(4)/(K(4,5)+R(4)),r(5)*R(5)/(K(5,5)+R(5))]);

N1dot = @(N, R) N(1)*(mu1(R) - m(1));
N2dot = @(N, R) N(2)*(mu2(R) - m(2));
N3dot = @(N, R) N(3)*(mu3(R) - m(3));
N4dot = @(N, R) N(4)*(mu4(R) - m(4));
N5dot = @(N, R) N(5)*(mu5(R) - m(5));

R1dot = @(N, R) D*(S(1)-R(1)) - sum([C(1,1)*mu1(R)*N(1),C(1,2)*mu2(R)*N(2), C(1,3)*mu3(R)*N(3), C(1,4)*mu4(R)*N(4),C(1,5)*mu5(R)*N(5)]);
R2dot = @(N, R) D*(S(2)-R(2)) - sum([C(2,1)*mu1(R)*N(2),C(2,2)*mu2(R)*N(2), C(2,3)*mu3(R)*N(3), C(2,4)*mu4(R)*N(4),C(2,5)*mu5(R)*N(5)]);
R3dot = @(N, R) D*(S(3)-R(3)) - sum([C(3,1)*mu1(R)*N(2),C(3,2)*mu2(R)*N(2), C(3,3)*mu3(R)*N(3), C(3,4)*mu4(R)*N(4),C(3,5)*mu5(R)*N(5)]);
R4dot = @(N, R) D*(S(4)-R(4)) - sum([C(4,1)*mu1(R)*N(2),C(4,2)*mu2(R)*N(2), C(4,3)*mu3(R)*N(3), C(4,4)*mu4(R)*N(4),C(4,5)*mu5(R)*N(5)]);
R5dot = @(N, R) D*(S(5)-R(5)) - sum([C(5,1)*mu1(R)*N(2),C(5,2)*mu2(R)*N(2), C(5,3)*mu3(R)*N(3), C(5,4)*mu4(R)*N(4),C(5,5)*mu5(R)*N(5)]);

%% initial conditions next up: 

M = ceil(T/deltat)+1; % total number of iterates we consider

R = zeros(M, size(K,1)); % an array to store the resource availabilities

N = zeros(M, size(K,2)); % an array to store the population abudances

time_vals = 0:deltat:(M-1)*deltat;

R0 = S; % initial conditions for the resources

N0 = (1/10) + (1/100).*(1:size(K,2));

R(1,:) = R0;

N(1,:) = N0;

for ii = 2:M

    % programming the (multivariate) fourth-order Runge-Kutta algorithm
    
    R1 = deltat * [R1dot(N0, R0), R2dot(N0, R0), R3dot(N0, R0), R4dot(N0, R0), R5dot(N0, R0)]; % evaluating the derivative at time tn
    
    N1 = deltat * [N1dot(N0, R0), N2dot(N0, R0), N3dot(N0, R0), N4dot(N0, R0), N5dot(N0, R0)]; % evaluating the derivative at time tn
    
    R2 = deltat * [R1dot(N0+N1/2, R0+R1/2), R2dot(N0+N1/2, R0+R1/2), R3dot(N0+N1/2, R0+R1/2), R4dot(N0+N1/2, R0+R1/2), R5dot(N0+N1/2, R0+R1/2)]; % evaluating the derivative at time tn + deltat/2
    
    N2 = deltat * [N1dot(N0+N1/2, R0+R1/2), N2dot(N0+N1/2, R0+R1/2), N3dot(N0+N1/2, R0+R1/2), N4dot(N0+N1/2, R0+R1/2), N5dot(N0+N1/2, R0+R1/2)]; % evaluating the derivative at time tn + deltat/2

    R3 = deltat * [R1dot(N0+N2/2, R0+R2/2), R2dot(N0+N2/2, R0+R2/2), R3dot(N0+N2/2, R0+R2/2), R4dot(N0+N2/2, R0+R2/2), R5dot(N0+N2/2, R0+R2/2)]; % evaluating the derivative at time tn + deltat/2
        
    N3 = deltat * [N1dot(N0+N2/2, R0+R2/2), N2dot(N0+N2/2, R0+R2/2), N3dot(N0+N2/2, R0+R2/2), N4dot(N0+N2/2, R0+R2/2), N5dot(N0+N2/2, R0+R2/2)]; % evaluating the derivative at time tn + deltat/2
        
    R4 = deltat * [R1dot(N0+N3, R0+R3), R2dot(N0+N3, R0+R3), R3dot(N0+N3, R0+R3), R4dot(N0+N3, R0+R3), R5dot(N0+N3, R0+R3)]; % evaluating the derivative at time tn + deltat

    N4 = deltat * [N1dot(N0+N3, R0+R3), N2dot(N0+N3, R0+R3), N3dot(N0+N3, R0+R3), N4dot(N0+N3, R0+R3), N5dot(N0+N3, R0+R3)]; % evaluating the derivative at time tn + deltat
        
    Rtot = (R1 + 2*R2 + 2*R3 + R4)/6; % taking the average
       
    Ntot = (N1 + 2*N2 + 2*N3 + N4)/6; % taking the average

    R0 = R0 + Rtot; % updating the value of the x variable
        
    N0 = N0 + Ntot; % updating the value of the y variable
    
    R(ii,:) = R0; % storing the new values
    N(ii,:) = N0; % storing the new values
    
end

% for ii = 2:M
%     % forward euler
%     R0 = R0 + deltat*[R1dot(N0, R0), R2dot(N0, R0), R3dot(N0, R0), R4dot(N0, R0), R5dot(N0, R0)];
%     N0 = N0 + deltat*[N1dot(N0, R0), N2dot(N0, R0), N3dot(N0, R0), N4dot(N0, R0), N5dot(N0, R0)]; 
%     R(ii,:) = R0; % storing the new values
%     N(ii,:) = N0; % storing the new values
% end

figure(1)

plot(time_vals, N(:,1))
hold on
plot(time_vals, N(:,2))
plot(time_vals, N(:,3))
plot(time_vals, N(:,4))
plot(time_vals, N(:,5))
xlabel('Time, days')
ylabel('Species abundances')
legend('Species 1', 'Species 2', 'Species 3', 'Species 4', 'Species 5', 'location', 'southeast')

figure(2)

plot3(N(:,1), N(:,3), N(:,5), 'k', 'LineWidth', 1)
hold on
%plot3(N(floor(M/2):M,1), N(floor(M/2):M,2), N(floor(M/2):M,3), '--g', 'LineWidth', 1.5)
xlabel('Species 1')
ylabel('Species 3')
zlabel('Species 5')
%title('Three-dimensional phase portrait')
%legend('Example trajectory', 'Limit cycle')

figure(3)

plot(time_vals, sum(N,2))
xlabel('Time, days')
ylabel('Total biomass')
title('Total biomass against time')