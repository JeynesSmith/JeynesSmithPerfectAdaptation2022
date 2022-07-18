% The script is used to generate heatmap plots for Figure 4 of the manuscript.
% This script generates abundance data for varying Etot and Itot under the
% complex-complete framework. Network parameters are set in the first half
% of the script.

%% Network Parameters

% Set the initial values
Atot = 200; Btot = 200; % total substrate abundances
C1 = 0; C2 = 0; C3 = 0; C4 = 0;
A = Atot; As = Atot-A;
B = Btot; Bs = Btot-B;
init = [A As B Bs C1 C2 C3 C4]; % initial concentration vector

% Set the model parameters
k1 = 200; k2 = 200; k3 = 10; k4 = 4; % catalytic constants
K1 = 1e0; K2 = 1e0; K3 = 1e-2; K4 = 1e-2; % Michaelis constants
d1=1; d2=1; d3=1; d4=1;
a1=(k1+d1)/K1;a2=(k2+d2)/K2;a3=(k3+d3)/K3;a4=(k4+d4)/K4;
k = [a1 d1 k1 a2 d2 k2 a3 d3 k3 a4 d4 k4]; % parameter vector

% Set network inputs
ET = 10.^[-3:3]; % total enzyme vector
IT = 10.^[-7:3]; % total input vector

%% simulation
tspan = [0 1e6]; % time span
Abundance = zeros(length(ET),length(IT),4); % initialise output abundance storage
sigma = k4./k3.*ET; % estimated set point

for i = 1:length(ET)
    for j = 1:length(IT)
        
        initAll = [init(1:4) IT(j) ET(i) init(5:end)]; % update initial condition for simulation
        
        % Perform the numerical integration for closed loop
        [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, initAll);
        
        % store relative abundances
        Abundance(i,j,1) = log10(abs(sigma(i)-sum(u(end,[2,9])))./sigma(i));
        Abundance(i,j,2) = log10(abs(sigma(i)-u(end,2))./sigma(i));
        Abundance(i,j,3) = log10(abs(sigma(i)-u(end,9))./sigma(i));
        Abundance(i,j,4) = log10(abs(ET(i)-u(end,10))./ET(i));
        
    end
end

% plot sigma vs output (A*_S)
figure
imagesc(Abundance(:,:,1))
caxis([-4,3])
a = colorbar;
a.Label.String = 'log_{10}(|\sigma - A^*_S|/\sigma)';
axis equal
xlabel('log_{10}(I_{tot})')
ylabel('log_{10}(E_{tot})')

% plot sigma vs A*
figure
imagesc(Abundance(:,:,2))
caxis([-4,3])
a = colorbar;
a.Label.String = 'log_{10}(|\sigma - A^*|/\sigma)';
axis equal
xlabel('log_{10}(I_{tot})')
ylabel('log_{10}(E_{tot})')

% plot sigma vs C3
figure
imagesc(Abundance(:,:,3))
caxis([-4,3])
a = colorbar;
a.Label.String = 'log_{10}(|\sigma - C_3|/\sigma)';
axis equal
xlabel('log_{10}(I_{tot})')
ylabel('log_{10}(E_{tot})')

% plot Etot vs C4
figure
imagesc(Abundance(:,:,4))
caxis([-4,3])
a = colorbar;
a.Label.String = 'log_{10}(|E_{tot} - C_4|/E_{tot})';
axis equal
xlabel('log_{10}(I_{tot})')
ylabel('log_{10}(E_{tot})')


%% ode system
function eqns = odesys(t,u,k)
eqns = zeros(10,1); % To start with we have twelve empty equations
% Using u = [C  Cs  B Bs  I E1 C1 C2 C3 C4]
% Using k = [a1,d1,k1,a2,d2,k2,a3,d3,k3,a4,d4,k4]
eqns(1) = k(2)*u(7) + k(6)*u(8) - k(1)*u(1)*u(5);
eqns(2) = k(3)*u(7) + k(5)*u(8) + k(8)*u(9) + k(9)*u(9) - k(7)*u(2)*u(3)  - k(4)*u(2)*u(4);
eqns(3) = k(8)*u(9) + k(12)*u(10) - k(7)*u(3)*u(2);
eqns(4) = k(5)*u(8) + k(6)*u(8) + k(11)*u(10) + k(9)*u(9) - k(4)*u(2)*u(4) - k(10)*u(4)*u(6);
eqns(5) = (k(2) + k(3))*u(7) - k(1)*u(1)*u(5);
eqns(6) = (k(11) + k(12))*u(10) - k(10)*u(4)*u(6);
eqns(7) = k(1)*u(1)*u(5) - (k(2) + k(3))*u(7);
eqns(8) = k(4)*u(2)*u(4) - (k(5) + k(6))*u(8);
eqns(9) = k(7)*u(2)*u(3) - (k(8) + k(9))*u(9);
eqns(10) = k(10)*u(4)*u(6) - (k(11) + k(12))*u(10);
end