% Set the initial values
Wtot = 100;
C1 = 0; C2 = 0;
W = Wtot; Ws = 0;
E2 = 20;

% Set the model parameters
k1 = 1; k2 = 1;
K2 = 100; K1 = K2;
d1=1;d2=1;a1=(1+k1)/K1;a2=(1+k2)/K2;
k = [a1,d1,k1,a2,d2,k2];

% Start Looping Process
x= 0:0.1:40; % set inputs
Wsplot = zeros(size(x));
for i = 1:length(x)
    
    E1 = x(i);
    init = [W Ws E1 E2 C1 C2];
    
    % Set the time domain
    tspan = [0 1e6];
    
    % Perform the numerical integration
    [t,u] = ode23tb(@(t,u) odesys(t,u,k), tspan, init);
    
    % Create vectors for plots
    Wsplot(i) = u(end,2)/Wtot;
end

% Plots
figure(1), hold on
plot(x,Wsplot,'-','LineWidth',2)
ylim([0,1])


%%ODE System
function eqns = odesys(t,u,k)
eqns = zeros(6,1); % To start with we have six empty equations
% Using u = [W  Ws E1 E2 C1 C2] 
% Using k = [a1,d1,k1,a2,d2,k2]
%           [ 1  2  3  4  5  6 
% k([3,6]) = k([3,6])+0.001.*(rand(1,2)-0.5);
eqns(1) = k(2)*u(5) + k(6)*u(6) - k(1)*u(1)*u(3);
eqns(2) = k(5)*u(6) + k(3)*u(5) - k(4)*u(2)*u(4);
eqns(3) = k(2)*u(5) + k(3)*u(5) - k(1)*u(1)*u(3);
eqns(4) = (k(5) + k(6))*u(6) - k(4)*u(2)*u(4);
eqns(5) = k(1)*u(1)*u(3) - (k(2) + k(3))*u(5);
eqns(6) = k(4)*u(2)*u(4) - (k(5) + k(6))*u(6);
end
