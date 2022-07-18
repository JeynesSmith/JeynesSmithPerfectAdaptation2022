% The script is used to generate simulations of a single reverserible
% covalent modification cycle for Figures 9 (altering E2) and 13 (altering 
% catalytic and Michaelis constants) of the manuscript. This script 
% generates abundance data for varying input, x. Network parameters are 
% set in the first half of the script.

%% Network Parameters

% Set the initial values
Wtot = 100; % total substrate
C1 = 0; C2 = 0;
W = Wtot; Ws = 0;
E2 = 20; % second enzyme abundance (B*_S)
x= 0:0.1:40; % set inputs

% Set the model parameters
k1 = 1; k2 = 1; % catalytic constants
K2 = 100; K1 = K2; % Michaelis constants
d1=1;d2=1;a1=(1+k1)/K1;a2=(1+k2)/K2;
k = [a1,d1,k1,a2,d2,k2];

%% Start Looping Process
Wsplot = zeros(size(x)); % initialise output vector

for i = 1:length(x)
    
    init = [W Ws x(i) E2 C1 C2];
    
    % Set the time domain
    tspan = [0 1e6];
    
    % Perform the numerical integration
    [t,u] = ode23tb(@(t,u) odesys(t,u,k), tspan, init);
    
    % Create vectors for plots
    Wsplot(i) = u(end,2)/Wtot; % abundance relative to total
end

% Plots
figure(1), hold on % plots on top of Figure 1 without clearing
plot(x,Wsplot,'-','LineWidth',2)
ylim([0,1])


%% ODE System
function eqns = odesys(t,u,k)
eqns = zeros(6,1); % To start with we have six empty equations
% Using u = [W  Ws E1 E2 C1 C2] 
% Using k = [a1,d1,k1,a2,d2,k2]
%           [ 1  2  3  4  5  6 
eqns(1) = k(2)*u(5) + k(6)*u(6) - k(1)*u(1)*u(3);
eqns(2) = k(5)*u(6) + k(3)*u(5) - k(4)*u(2)*u(4);
eqns(3) = k(2)*u(5) + k(3)*u(5) - k(1)*u(1)*u(3);
eqns(4) = (k(5) + k(6))*u(6) - k(4)*u(2)*u(4);
eqns(5) = k(1)*u(1)*u(3) - (k(2) + k(3))*u(5);
eqns(6) = k(4)*u(2)*u(4) - (k(5) + k(6))*u(6);
end
