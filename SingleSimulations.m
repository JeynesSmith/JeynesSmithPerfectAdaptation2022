% The script is used to generate single simulations of the network under
% the multiple different frameworks: Ma et al. and Ferrell for Figure 3
% (MaPlot,FerrellPlot), Ma et al. vs complex-complete for Figure 5
% (MaCCPlot), and the complex-complete model for Figure 6 (CCSimplePlot).
% This script generates the data used in the figures, with the parameters
% being set in the first section.

%% Parameters and model selection
% Set the initial values
Atot = 1e3; Btot = 1e4;
Etot = 10^1.5;
C1 = 0; C2 = 0; C3 = 0; C4 = 0;
C = Atot; Cs = Atot-C;
B = Btot; Bs = Btot-B;
init = [C Cs B Bs Etot C1 C2 C3 C4];

% Set the model parameters
k1 = 200; k2 = 200; k3 = 10; k4 = 4;
K1 = 1e0; K2 = 1e0; K3 = 1e-2; K4 = 1e-2;
d1=1; d2=1; d3=1; d4=1;
a1=(k1+d1)/K1;a2=(k2+d2)/K2;a3=(k3+d3)/K3;a4=(k4+d4)/K4;
k = [a1 d1 k1 a2 d2 k2 a3 d3 k3 a4 d4 k4];

% Choose model to simulate
RangeCC = 0; % Calculate the range of the complexcomplete model
CCPlot = 1; % plot the complex-complete model
CCSimplePlot = 0; % plot the complex-complete model with simple outputs
RangeMa = 0; % Calculate the range of the Ma et al. model
MaPlot = 1; % plot the Ma et al. model
MaCCPlot = 1; % plot the Ma et al. model and complex-complete simple outputs
RangeFerrell = 0; % Calculate the range of the Ferrell model
FerrellPlot = 1; % plot the Ferrell model

% set simulation parameters
Atol = 1e-2; % set A change tolerance
Btol = 1e-2; % set B change tolerance
tspan = [0 1e6]; % time span
maxI = 1e15; % max accepted input
maxSim = 30; % max orders of magnitude after I_S
maxCount = 10; % max simulations for refined search
startI = 1e-5; % Initial input to start search from
inputStep = startI; % define step size
theoreticalEndPoint = k4/k3*Etot; % define the estimated setpoint
changeInput = [1 7/8]; % small input perturbation size

% formatting
format shortg, warning off

%% Complex-complete range finder
if RangeCC == 1 || CCSimplePlot == 1 || CCPlot == 1
    
    
    tempAbund = zeros(2,10); % initialise abundance matrix
    % broad search (increasing by order of magnitude at a time)
    for j = 1:2
        initAll = [init(1:4) startI/changeInput(j) init(5:end)];
        [t,u] = ode23s(@(t,u) odesys(t,u,k), tspan, initAll);
        tempAbund(j,:) = u(end,:);
    end
    
    while ~(abs(sum(tempAbund(1,[2,9])) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ... % check that A*S is close to setpoint
            && abs(sum(tempAbund(1,[2,9])) - sum(tempAbund(2,[2,9])))/abs(sum(tempAbund(1,[2,9])))<=2*Atol ... % check that A*S isn't changing
            && abs(sum(tempAbund(1,[4,8])) - sum(tempAbund(2,[4,8])))/abs(sum(tempAbund(1,[4,8])))>=Btol) ... % check that B*S is changing
            && maxI > startI % check that we haven't exceeded the maximum input to avoid infinite search
        startI = startI*10; % increase input by order of magnitude
        tempAbund(2,:) = tempAbund(1,:);
        initAll = [init(1:4) startI init(5:end)];
        [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, initAll);
        tempAbund(1,:) = u(end,:); % update abundance
    end
    % refined search (increasing by order of 0.1*startI from broad search)
    startI = startI/5; % decrease input by one order of magnitude and then increase
    stepSize = startI/2; % set step size
    count = 2; % initialise count
    for j = 1:2
        initAll = [init(1:4) startI/j init(5:end)];
        [t,u] = ode23s(@(t,u) odesys(t,u,k), tspan, initAll);
        tempAbund(j,:) = u(end,:);
    end
    % repeat checks and updates for smaller input increase
    while ~(abs(sum(tempAbund(1,[2,9])) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
            && abs(sum(tempAbund(1,[2,9])) - sum(tempAbund(2,[2,9])))/abs(sum(tempAbund(1,[2,9])))<=2*Atol ...
            && abs(sum(tempAbund(1,[4,8])) - sum(tempAbund(2,[4,8])))/abs(sum(tempAbund(1,[4,8])))>=Btol) ...
            && maxI > startI && count<maxCount % added check to see if returned to previous order of magnitude
        startI = startI + stepSize; % increase by stepsize
        tempAbund(2,:) = tempAbund(1,:);
        initAll = [init(1:4) startI init(5:end)];
        [t,u] = ode23s(@(t,u) odesys(t,u,k), tspan, initAll);
        tempAbund(1,:) = u(end,:);
        count = count + 1;
    end
    simNum = 0;
    finalI = startI; % set the final input to be the start input
    % repeat updates until we reach an order of magnitude or until RPA is lost
    while abs(sum(tempAbund(1,[2,9])) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
            && abs(sum(tempAbund(1,[2,9])) - sum(tempAbund(2,[2,9])))/abs(sum(tempAbund(1,[2,9])))<=2*Atol ...
            && abs(sum(tempAbund(1,[4,8])) - sum(tempAbund(2,[4,8])))/abs(sum(tempAbund(1,[4,8])))>Btol ...
            && maxI > finalI && count<maxCount
        finalI = finalI + stepSize;
        tempAbund(2,:) = tempAbund(1,:);
        initAll = [init(1:4) finalI init(5:end)];
        [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, initAll);
        tempAbund(1,:) = u(end,:);
        count = count + 1;
    end
    % check that we still have RPA
    if abs(sum(tempAbund(1,[2,9])) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
            && abs(sum(tempAbund(1,[2,9])) - sum(tempAbund(2,[2,9])))/abs(sum(tempAbund(1,[2,9])))<=2*Atol ...
            && abs(sum(tempAbund(1,[4,8])) - sum(tempAbund(2,[4,8])))/abs(sum(tempAbund(1,[4,8])))>Btol
        
        tempAbund(2,:) = tempAbund(1,:);
        finalI = finalI*10;
        initAll = [init(1:4) finalI init(5:end)];
        [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, initAll);
        tempAbund(1,:) = u(end,:);
        simNum = 1;
        % broad search (increasing by order of magnitude at a time)
            % repeat checks and updates
        while abs(sum(tempAbund(1,[2,9])) - theoreticalEndPoint)/abs(theoreticalEndPoint) <= Atol ...
                && abs(sum(tempAbund(1,[2,9])) - sum(tempAbund(2,[2,9])))/abs(sum(tempAbund(1,[2,9])))<=2*Atol ...
                && abs(sum(tempAbund(1,[4,8])) - sum(tempAbund(2,[4,8])))/abs(sum(tempAbund(1,[4,8]))) > Btol ...
                && maxSim >= simNum
            finalI = finalI*10;
            tempAbund(2,:) = tempAbund(1,:);
            initAll = [init(1:4) finalI init(5:end)];
            [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, initAll);
            tempAbund(1,:) = u(end,:);
            simNum = simNum + 1;
        end
        % refined search (increasing by order of 0.1*startI from broad search)
        stepSize = finalI/10; % set step size
        finalI = max(startI,finalI/5); % decrease input by one order of magnitude and then increase
        count = 2;
        for j = 1:2
            initAll = [init(1:4) finalI/j init(5:end)];
            [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, initAll);
            tempAbund(j,:) = u(end,:);
        end
        % repeat checks and updates for smaller stepsize
        while abs(sum(tempAbund(1,[2,9])) - theoreticalEndPoint)/abs(theoreticalEndPoint) <= Atol ...
                && abs(sum(tempAbund(1,[2,9])) - sum(tempAbund(2,[2,9])))/abs(sum(tempAbund(1,[2,9])))<=2*Atol ...
                && abs(sum(tempAbund(1,[4,8])) - sum(tempAbund(2,[4,8])))/abs(sum(tempAbund(1,[4,8]))) > Btol ...
                && count<maxCount
            finalI = finalI + stepSize;
            tempAbund(2,:) = tempAbund(1,:);
            initAll = [init(1:4) finalI init(5:end)];
            [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, initAll);
            tempAbund(1,:) = u(end,:);
            simNum = simNum + 1;
            count = count + 1;
        end
    end
    disp([startI, finalI])
end

%% Complex-complete plot
if CCPlot == 1 || CCSimplePlot == 1
    
    x = linspace(0,1.5*finalI,101); % define input vector across RPA Range
    % x = linspace(0,5,200); % manually set input vector
    
    Abundance = zeros(length(x),10); % initialise abundance matrix
    
    for i = 1:length(x) % loop over inputs
        
        initAll = [init(1:4) x(i) init(5:end)];
        
        % Perform the numerical integration
        [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan.^2, initAll);
        
        % store data for plot
        Abundance(i,:) = u(end,:);
        
    end
    
    % plot
    figure
    yyaxis left % RPA variable
    plot(x,sum(Abundance(:,[2,9]),2),'-','LineWidth',2,'color',[0,0.5,1]) % total input to opposer cycle
    hold on
    if CCSimplePlot == 1
        plot(x,Abundance(:,2),'--','LineWidth',2,'color',[0,0.4,0.9]) % A*
    end
    plot(x,ones(size(x)).*theoreticalEndPoint,'k:')% set point
    ylabel('Protein Concentration')
    yyaxis right % nonRPA variable
    plot(x,sum(Abundance(:,[4,8]),2),'-','LineWidth',2,'color',[1,0.5,0]) % total input to input/output cycle
    hold on
    if CCSimplePlot == 1
        plot(x,Abundance(:,4),'--','LineWidth',2,'color',[0.9,0.4,0]) % B*
    end
    ylabel('Protein Concentration')
    xlabel('Input (I_{tot})')
    
    if CCSimplePlot == 1
        legend('{A^*}_S','A^*','\sigma','{B^*}_S','B^*')
    else
        legend('{A^*}_S','\sigma','{B^*}_S')
    end
    
end


%% Ma et al. range finder
if RangeMa == 1 || MaPlot == 1 || MaCCPlot
    
    % define k (parameters) and init (initial conditions)
    init = [Atot Btot];
    MaxTot = max(init);
    kinit = [k1 k2 k3 k4 K1 K2 K3 K4 Etot Atot Btot];
    theoreticalEndPoint = kinit(4)/kinit(3)*kinit(9); % define the estimated setpoint
    
    % calculate first input for which RPA can be observed
    startIMa = 1e-5; % initialise input
    tempAbund = zeros(2,2); % initialise abundance matrix
    %broad search (increasing by order of magnitude at a time)
    for j = 1:2
        k = [kinit(1:8) startIMa/changeInput(j) kinit(9:end)];
        [~,u] = ode23s(@(t,u) odesysMa(t,u,k,MaxTot), tspan, init);
        tempAbund(j,:) = u(end,:);
    end
    while ~(abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ... % check if A* within tolerance of setpoint
            && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ... % check if A* changed within tolerance
            && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2))>=Btol) ... % check if B* changed more than tolerance
            && maxI > startIMa % stop if input gets too large
        startIMa = startIMa*10; % increase input by order of magnitude
        tempAbund(2,:) = tempAbund(1,:);
        k = [kinit(1:8) startIMa kinit(9:end)];
        [~,u] = ode23s(@(t,u) odesysMa(t,u,k,MaxTot), tspan, init);
        tempAbund(1,:) = u(end,:); % update abundance
    end
    % refined search  (increasing by order of 0.1*startI from broad search)
    startIMa = startIMa/5; % decrease input by one order of magnitude and then increase
    stepSize = startIMa/2; % set step size
    count = 2; % initialise count
    for j = 1:2
        k = [kinit(1:8) startIMa/j kinit(9:end)];
        [t,u] = ode23s(@(t,u) odesysMa(t,u,k,MaxTot), tspan, init);
        tempAbund(j,:) = u(end,:);
    end
    % repeat checks and updates for smaller input increase
    while ~(abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
            && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
            && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2))>=Btol) ...
            && maxI > startIMa && count<maxCount % added check to see if returned to previous order of magnitude
        startIMa = startIMa + stepSize; % increase by stepsize
        tempAbund(2,:) = tempAbund(1,:);
        k = [kinit(1:8) startIMa kinit(9:end)];
        [t,u] = ode23s(@(t,u) odesysMa(t,u,k,MaxTot), tspan, init);
        tempAbund(1,:) = u(end,:);
        count = count + 1;
    end
    simNum = 0;
    finalIMa = startIMa; % set the final input to be the start input
    % repeat updates until we reach an order of magnitude or until RPA is lost
    while abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
            && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
            && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2))>Btol ...
            && maxI > finalIMa && count<maxCount
        finalIMa = finalIMa + stepSize;
        tempAbund(2,:) = tempAbund(1,:);
        k = [kinit(1:8) finalIMa kinit(9:end)];
        [~,u] = ode23s(@(t,u) odesysMa(t,u,k,MaxTot), tspan, init);
        tempAbund(1,:) = u(end,:);
        count = count + 1;
    end
    
    % check that we still have RPA
    if abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
            && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
            && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2))>Btol
        
        tempAbund(2,:) = tempAbund(1,:);
        finalIMa = finalIMa*10;
        k = [kinit(1:8) finalIMa kinit(9:end)];
        [~,u] = ode23s(@(t,u) odesysMa(t,u,k,MaxTot), tspan, init);
        tempAbund(1,:) = u(end,:);
        simNum = 1;
        % broad search (increasing by order of magnitude at a time)
        % repeat checks and updates
        while abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint) <= Atol ...
                && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
                && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2)) > Btol ...
                && maxSim >= simNum
            finalIMa = finalIMa*10;
            tempAbund(2,:) = tempAbund(1,:);
            k = [kinit(1:8) finalIMa kinit(9:end)];
            [~,u] = ode23s(@(t,u) odesysMa(t,u,k,MaxTot), tspan, init);
            tempAbund(1,:) = u(end,:);
            simNum = simNum + 1;
        end
        % refined search (increasing by order of 0.1*startI from broad search)
        stepSize = finalIMa/10; % set step size
        finalIMa = max(startIMa,finalIMa/5); % decrease input by one order of magnitude and then increase
        count = 2;
        for j = 1:2
            k = [kinit(1:8) finalIMa/j kinit(9:end)];
            [~,u] = ode23s(@(t,u) odesysMa(t,u,k,MaxTot), tspan, init);
            tempAbund(j,:) = u(end,:);
        end
        % repeat checks and updates for smaller stepsize
        while abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint) <= Atol ...
                && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
                && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2)) > Btol ...
                && count<maxCount
            finalIMa = finalIMa + stepSize;
            tempAbund(2,:) = tempAbund(1,:);
            k = [kinit(1:8) finalIMa kinit(9:end)];
            [~,u] = ode23s(@(t,u) odesysMa(t,u,k,MaxTot), tspan, init);
            tempAbund(1,:) = u(end,:);
            simNum = simNum + 1;
            count = count + 1;
        end
    end
    disp([startIMa,finalIMa])
end

%% Ferrell range finder
if RangeFerrell == 1 || FerrellPlot
    
    % define k (parameters) and init (initial conditions)
    init = [Atot Btot];
    MaxTot = max(init);
    kinit = [k1 k2 k3 k4 K1 K2 K3 K4 Etot Atot Btot];
    theoreticalEndPoint = kinit(4)/kinit(3)*kinit(9); % define the estimated setpoint
    
    % calculate first input for which RPA can be observed
    startIFerrell = 1e-8;
    tempAbund = zeros(2,2);
    %broad search (increasing by order of magnitude at a time)
    for j = 1:2
        k = [kinit(1:8) startIFerrell/changeInput(j) kinit(9:end)];
        [~,u] = ode23s(@(t,u) odesysFerrell(t,u,k,MaxTot), tspan, init);
        tempAbund(j,:) = u(end,:);
    end
    while ~(abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ... % check if A* within tolerance of setpoint
            && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ... % check if A* changed within tolerance
            && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2))>=Btol) ... % check if B* changed more than tolerance
            && maxI > startIFerrell % stop if input gets too large
        startIFerrell = startIFerrell*10; % increase input by order of magnitude
        tempAbund(2,:) = tempAbund(1,:);
        k = [kinit(1:8) startIFerrell kinit(9:end)];
        [~,u] = ode23s(@(t,u) odesysFerrell(t,u,k,MaxTot), tspan, init);
        tempAbund(1,:) = u(end,:); % update abundance
    end
    % refined search (increasing by order of 0.1*startI from broad search)
    startIFerrell = startIFerrell/5; % decrease input by one order of magnitude and then increase
    stepSize = startIFerrell/2; % set step size
    count = 2; % initialise count
    for j = 1:2
        k = [kinit(1:8) startIFerrell/j kinit(9:end)];
        [t,u] = ode23s(@(t,u) odesysFerrell(t,u,k,MaxTot), tspan, init);
        tempAbund(j,:) = u(end,:);
    end
    % repeat checks and updates for smaller input increase
    while ~(abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
            && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
            && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2))>=Btol) ...
            && maxI > startIFerrell && count<maxCount % added check to see if returned to previous order of magnitude
        startIFerrell = startIFerrell + stepSize; % increase by stepsize
        tempAbund(2,:) = tempAbund(1,:);
        k = [kinit(1:8) startIFerrell kinit(9:end)];
        [t,u] = ode23s(@(t,u) odesysFerrell(t,u,k,MaxTot), tspan, init);
        tempAbund(1,:) = u(end,:);
        count = count + 1;
    end
    simNum = 0;
    finalIFerrell = startIFerrell; % set the final input to be the start input
    % repeat updates until we reach an order of magnitude or until RPA is lost
    while abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
            && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
            && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2))>Btol ...
            && maxI > finalIFerrell && count<maxCount
        finalIFerrell = finalIFerrell + stepSize;
        tempAbund(2,:) = tempAbund(1,:);
        k = [kinit(1:8) finalIFerrell kinit(9:end)];
        [~,u] = ode23s(@(t,u) odesysFerrell(t,u,k,MaxTot), tspan, init);
        tempAbund(1,:) = u(end,:);
        count = count + 1;
    end
    
    % check that we still have RPA
    if abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
            && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
            && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2))>Btol
        
        tempAbund(2,:) = tempAbund(1,:);
        finalIFerrell = finalIFerrell*10;
        k = [kinit(1:8) finalIFerrell kinit(9:end)];
        [~,u] = ode23s(@(t,u) odesysFerrell(t,u,k,MaxTot), tspan, init);
        tempAbund(1,:) = u(end,:);
        simNum = 1;
        % broad search (increasing by order of magnitude at a time)
        % repeat checks and updates
        while abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint) <= Atol ...
                && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
                && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2)) > Btol ...
                && maxSim >= simNum
            finalIFerrell = finalIFerrell*10;
            tempAbund(2,:) = tempAbund(1,:);
            k = [kinit(1:8) finalIFerrell kinit(9:end)];
            [~,u] = ode23s(@(t,u) odesysFerrell(t,u,k,MaxTot), tspan, init);
            tempAbund(1,:) = u(end,:);
            simNum = simNum + 1;
        end
        % refined search (increasing by order of 0.1*startI from broad search)
        stepSize = finalIFerrell/10; % set step size
        finalIFerrell = max(startIFerrell,finalIFerrell/5); % decrease input by one order of magnitude and then increase
        count = 2;
        for j = 1:2
            k = [kinit(1:8) finalIFerrell/j kinit(9:end)];
            [~,u] = ode23s(@(t,u) odesysFerrell(t,u,k,MaxTot), tspan, init);
            tempAbund(j,:) = u(end,:);
        end
        % repeat checks and updates for smaller stepsize
        while abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint) <= Atol ...
                && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
                && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2)) > Btol ...
                && count<maxCount
            finalIFerrell = finalIFerrell + stepSize;
            tempAbund(2,:) = tempAbund(1,:);
            k = [kinit(1:8) finalIFerrell kinit(9:end)];
            [~,u] = ode23s(@(t,u) odesysFerrell(t,u,k,MaxTot), tspan, init);
            tempAbund(1,:) = u(end,:);
            simNum = simNum + 1;
            count = count + 1;
        end
    end
    disp([startIFerrell,finalIFerrell])
end

%% Plot Ma et al. Model
if MaCCPlot == 1 || MaPlot == 1
    
    x = linspace(0,1.5*finalIMa,100); % set input vector across RPA Range
    % x = linspace(0,5,200); % manually set input vector
    AbundanceMa = zeros(length(x),2); % initialise abundance matrix (Ma)
    Abundance = zeros(length(x),10); % initialise abundance matrix (CC)
    
    init = [Atot 0 Btot 0 Etot 0 0 0 0]; % initial conditions CC
    k = [a1 d1 k1 a2 d2 k2 a3 d3 k3 a4 d4 k4];
    kinit = [k1 k2 k3 k4 K1 K2 K3 K4 Etot Atot Btot]; % initial conditions Ma
    initSimp = [Atot Btot];
    
    % simulate
    for i = 1:length(x) % loop over input
        % Perform the numerical integration
        % CC model
        [~,u] = ode23s(@(t,u) odesys(t,u,k), [0 1e7], [init(1:4) x(i) init(5:end)]);
        Abundance(i,:) = u(end,:); % store CC abundances
        % Ma et al. model
        [~,u] = ode23s(@(t,u) odesysMa(t,u,[kinit(1:8) x(i) kinit(9:end)]), [0 1e7], initSimp);
        AbundanceMa(i,:) = u(end,:); % store Ma abundances
        
    end
    
    % plot
    figure
    if MaCCPlot == 1 % plot both CC and Ma
        yyaxis left % RPA variable
        plot(x,AbundanceMa(:,1),'-','LineWidth',2,'color',[0,0.5,1]) % ma output, A*
        hold on
        plot(x,Abundance(:,2),'--','LineWidth',2,'color',[0,0.4,0.9]) % CC output (A*)
        plot(x,ones(size(x)).*theoreticalEndPoint,'k:') % set point
        ylabel('Protein Concentration')
        yyaxis right % nonRPA variable
        plot(x,AbundanceMa(:,2),'-','LineWidth',2,'color',[1,0.5,0]) % ma output, B*
        hold on
        plot(x,Abundance(:,4),'--','LineWidth',2,'color',[0.9,0.4,0]) % CC output (B*)
        ylabel('Protein Concentration')
        xlabel('Input (I_{tot})')
        legend('Ma: A^*','CC: A^*','\sigma','Ma: B^*','CC: B^*')
    else % just plot Ma
        yyaxis left
        plot(x,AbundanceMa(:,1),'-','LineWidth',2,'color',[0,0.5,1]) % A*
        hold on
        plot(x,ones(size(x)).*theoreticalEndPoint,'k:') % set point
        ylabel('Protein Concentration')
        yyaxis right
        plot(x,AbundanceMa(:,2),'-','LineWidth',2,'color',[1,0.5,0]) % B*
        ylabel('Protein Concentration')
        xlabel('Input (I_{tot})')
        legend('A^*','\sigma','B^*')
    end
    
end

%% Plot Ferrell model
if FerrellPlot == 1
    
    x = linspace(0,1.5*finalIFerrell,101);  % set input vector across RPA Range
    % x = linspace(0,5,200); % manually set input vector
    AbundanceFerrell = zeros(length(x),2); % initialise abundance matrix
    
    kinit = [k1 k2 k3 k4 K1 K2 K3 K4 Etot Atot Btot]; % initial conditions
    initSimp = [Atot Btot];MaxTot = max(initSimp);
    
    % simulate
    for i = 1:length(x) % loop over input
        % Perform the numerical integration
        [~,u] = ode23s(@(t,u) odesysFerrell(t,u,[kinit(1:8) x(i) kinit(9:end)]), [0 1e7], initSimp);
        AbundanceFerrell(i,:) = u(end,:); % store Ferrell abundances
        
    end
    
    % plot
    figure
    yyaxis left
    plot(x,AbundanceFerrell(:,1),'-','LineWidth',2,'color',[0,0.5,1]) % A*
    hold on
    plot(x,ones(size(x)).*theoreticalEndPoint,'k:') % setpoint
    ylabel('Protein Concentration')
    yyaxis right
    plot(x,AbundanceFerrell(:,2),'-','LineWidth',2,'color',[1,0.5,0]) % B*
    ylabel('Protein Concentration')
    xlabel('Input (I_{tot})')
    legend('A^*','\sigma','B^*')
    
end

%% Complex-complete ODE system
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
%% Ferrell ODE System
function eqns = odesysFerrell(t,u,k,MaxTot)
eqns = zeros(2,1); % To start with we have twelve empty equations
% Using u = [A B]
% Using k = [k1,k2,k3,k4,K1,K2,K3,K4,I,E1,At,Bt]
u(u<0) = 0;
% u(u>1.1*MaxTot) = MaxTot;
eqns(1) = k(1)*k(9)*(k(11) - u(1)) - k(2)*u(1)*u(2);
eqns(2) = k(3)*u(1)*(k(12)-u(2))/(k(7) + k(12)-u(2)) - k(4)*k(10)*u(2)/(k(8) + u(2));
end
%% Ma et al. ODE System
function eqns = odesysMa(t,u,k,MaxTot)
eqns = zeros(2,1); % To start with we have twelve empty equations
% Using u = [A B]
% Using k = [k1,k2,k3,k4,K1,K2,K3,K4,I,E1,At,Bt]
u(u<0) = 0;
% u(u>1.1*MaxTot) = MaxTot;
eqns(1) = k(1)*k(9)*(k(11)-u(1))/(k(5) + k(11)-u(1)) - k(2)*u(1)*u(2)/(k(6) + u(1));
eqns(2) = k(3)*u(1)*(k(12)-u(2))/(k(7) + k(12)-u(2)) - k(4)*k(10)*u(2)/(k(8) + u(2));
end