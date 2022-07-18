% This script is used to generate data regarding how parameters affect the
% ability for the complex-complete network to generate the RPA response,
% and the range of inputs over which it can be observed. The parameters
% which are varied are set at the top of the first section, and the save
% file name is changed near the bottom of the script. 

%% generate random samples
rng(2) % initialise the random number generator for repeatability

% define the number of samples
numberSamples = 1e5;

% Define parameters from samples
%%% Totally random
% Parameters = [10.^(rand(numberSamples,8)*7-3) 10.^(rand(numberSamples,8)*4)];
%%% Random Michaelis constants
% Parameters = [ones(numberSamples,4).*[7,5,2,3] 10.^(rand(numberSamples,4)*7-3) ones(numberSamples,3).*[2e1 1e1 1e0]];
%%% Random catalytic constants
% Parameters = [ 10.^(rand(numberSamples,4)*7-3) ones(numberSamples,4).*[7e0,5e0,2e-2,3e-2] ones(numberSamples,3).*[2e1 1e1 1e0]];
%%% Random Total abundances
Parameters = [ ones(numberSamples,4).*[7,5,2,3] ones(numberSamples,4).*[7e0,5e0,2e-2,3e-2] 10.^(rand(numberSamples,3).*4)];

% simulation parameters
Atol = 1e-2; % set A change tolerance
Btol = 1e-2; % set B change tolerance
changeInput = [1 7/8]; % small input perturbation size
maxCount = 10; % max simulations for refined search
maxI = 1e7; % max accepted input
maxSim = 30; % max orders of magnitude after I_S
tspan = [0 1e6]; % max orders of magnitude after I_S

RPASets = zeros(numberSamples,size(Parameters,2)+2); % initialise RPA parameter sets
nonRPASets = zeros(numberSamples,size(Parameters,2)); % initialise nonRPA parameter sets

%% looping procedure for testing parameter sets
warning off
global EventOccur

for i = 1:numberSamples % loop over samples
    
    disp([num2str(i/numberSamples*100) '%'])
    
    % reset event check
    EventOccur = 0;
    
    % define k (parameters) and init (initial conditions)
    init = zeros(1,9);
    init([1,3,5]) = Parameters(i,9:11);
    MaxTot = max(Parameters(i,9:11)); % rough check for exceeding maxmimum protein abundance
    k = ones(1,12);
    k(3:3:12) = Parameters(i,1:4);
    k(1:3:10) = (k(3:3:12) + 1)./Parameters(i,5:8);
    theoreticalEndPoint = k(12)/k(9)*init(5); % define the estimated setpoint
    
    % calculate first input for which RPA can be observed
    startI = 1e-5; % initialise input
    tempAbund = zeros(2,10); % initialise abundance matrix
    % broad search (increasing by order of magnitude at a time)
    for j = 1:2
        initAll = [init(1:4) startI/changeInput(j) init(5:end)];
        tstart=tic;
        [~,u] = ode23s(@(t,u) odesys(t,u,k,max([MaxTot,initAll(5)])), tspan, initAll, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
        tempAbund(j,:) = u(end,:);
    end
    while ~(abs(sum(tempAbund(1,[2,9])) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ... % check that A*S is close to setpoint
            && abs(sum(tempAbund(1,[2,9])) - sum(tempAbund(2,[2,9])))/abs(sum(tempAbund(1,[2,9])))<=2*Atol ... % check that A*S isn't changing
            && abs(sum(tempAbund(1,[4,8])) - sum(tempAbund(2,[4,8])))/abs(sum(tempAbund(1,[4,8])))>=Btol) ...% check that B*S is changing
            && maxI > startI % check that we haven't exceeded the maximum input to avoid infinite search
        startI = startI*10; % increase input by order of magnitude
        tempAbund(2,:) = tempAbund(1,:); 
        initAll = [init(1:4) startI init(5:end)];
        tstart=tic; % initialise timer for time check
        [~,u] = ode23s(@(t,u) odesys(t,u,k,max([MaxTot,initAll(5)])), tspan, initAll, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
        tempAbund(1,:) = u(end,:); % update abundance
    end
    % refined search (increasing by order of 0.1*startI from broad search)
    startI = startI/5; % decrease input by one order of magnitude and then increase
    stepSize = startI/2; % set step size
    count = 2; % initialise count
    for j = 1:2
        initAll = [init(1:4) startI/j init(5:end)];
        tstart=tic;
        [t,u] = ode23s(@(t,u) odesys(t,u,k,max([MaxTot,initAll(5)])), tspan, initAll, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
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
        tstart=tic;
        [t,u] = ode23s(@(t,u) odesys(t,u,k,max([MaxTot,initAll(5)])), tspan, initAll, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
        tempAbund(1,:) = u(end,:);
        count = count + 1;
    end
    
    % if RPA, calculate end input
    if maxI ~= startI % check that a starting input was found
        
        finalInput = startI; % set the final input to be the start input
        % repeat updates until we reach an order of magnitude or until RPA is lost
        while abs(sum(tempAbund(1,[2,9])) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
                && abs(sum(tempAbund(1,[2,9])) - sum(tempAbund(2,[2,9])))/abs(sum(tempAbund(1,[2,9])))<=2*Atol ...
                && abs(sum(tempAbund(1,[4,8])) - sum(tempAbund(2,[4,8])))/abs(sum(tempAbund(1,[4,8])))>=Btol ...
                && maxI > finalInput && count<maxCount
            finalInput = finalInput + stepSize;
            tempAbund(2,:) = tempAbund(1,:);
            initAll = [init(1:4) finalInput init(5:end)];
            tstart=tic;
            [~,u] = ode23s(@(t,u) odesys(t,u,k,max([MaxTot,initAll(5)])), tspan, initAll, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
            tempAbund(1,:) = u(end,:);
            count = count + 1;
        end
        
        % check that we still have RPA
        if abs(sum(tempAbund(1,[2,9])) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
                && abs(sum(tempAbund(1,[2,9])) - sum(tempAbund(2,[2,9])))/abs(sum(tempAbund(1,[2,9])))<=2*Atol ...
                && abs(sum(tempAbund(1,[4,8])) - sum(tempAbund(2,[4,8])))/abs(sum(tempAbund(1,[4,8])))>=Btol
            
            tempAbund(2,:) = tempAbund(1,:);
            finalInput = finalInput*10;
            initAll = [init(1:4) finalInput init(5:end)];
            tstart=tic;
            [~,u] = ode23s(@(t,u) odesys(t,u,k,max([MaxTot,initAll(5)])), tspan, initAll, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
            tempAbund(1,:) = u(end,:);
            simNum = 1;
            % broad search (increasing by order of magnitude at a time)
            % repeat checks and updates
            while abs(sum(tempAbund(1,[2,9])) - theoreticalEndPoint)/abs(theoreticalEndPoint) <= Atol ...
                    && abs(sum(tempAbund(1,[2,9])) - sum(tempAbund(2,[2,9])))/abs(sum(tempAbund(1,[2,9])))<=2*Atol ...
                    && abs(sum(tempAbund(1,[4,8])) - sum(tempAbund(2,[4,8])))/abs(sum(tempAbund(1,[4,8]))) >= Btol ...
                    && maxSim >= simNum
                finalInput = finalInput*10;
                tempAbund(2,:) = tempAbund(1,:);
                initAll = [init(1:4) finalInput init(5:end)];
                tstart=tic;
                [~,u] = ode23s(@(t,u) odesys(t,u,k,max([MaxTot,initAll(5)])), tspan, initAll, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
                tempAbund(1,:) = u(end,:);
                simNum = simNum + 1;
            end
            % refined search (increasing by order of 0.1*startI from broad search)
            stepSize = finalInput/10; % set step size
            finalInput = max(startI,finalInput/5); % decrease input by one order of magnitude and then increase
            count = 2;
            for j = 1:2
                initAll = [init(1:4) finalInput/j init(5:end)];
                tstart=tic;
                [~,u] = ode23s(@(t,u) odesys(t,u,k,max([MaxTot,initAll(5)])), tspan, initAll, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
                tempAbund(j,:) = u(end,:);
            end
            % repeat checks and updates for smaller stepsize
            while abs(sum(tempAbund(1,[2,9])) - theoreticalEndPoint)/abs(theoreticalEndPoint) <= Atol ...
                    && abs(sum(tempAbund(1,[2,9])) - sum(tempAbund(2,[2,9])))/abs(sum(tempAbund(1,[2,9])))<=2*Atol ...
                    && abs(sum(tempAbund(1,[4,8])) - sum(tempAbund(2,[4,8])))/abs(sum(tempAbund(1,[4,8]))) >= Btol ...
                    && count<maxCount
                finalInput = finalInput + stepSize;
                tempAbund(2,:) = tempAbund(1,:);
                initAll = [init(1:4) finalInput init(5:end)];
                tstart=tic;
                [~,u] = ode23s(@(t,u) odesys(t,u,k,max([MaxTot,initAll(5)])), tspan, initAll, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
                tempAbund(1,:) = u(end,:);
                simNum = simNum + 1;
                count = count + 1;
            end
        end
        if EventOccur == 1 % check if simulation ever took too long to run
            disp('had to break')
            RPASets(i,:) = [Parameters(i,:) -1 -1]; % save in RPASets but with negative values of Is and If
        else
            % save rpa set parameters, first & final inputs
            disp('yes')
            RPASets(i,:) = [Parameters(i,:) startI finalInput];
        end
        
    elseif EventOccur == 1 % check if simulation ever took too long to run
        disp('had to break')
        RPASets(i,:) = [Parameters(i,:) -1 -1]; % save in RPASets but with negative values of Is and If
    else
        % save non rpa set parameters
        disp('no')
        nonRPASets(i,:) = Parameters(i,:);
    end
    
end
%% save data set
save CCRSTot nonRPASets RPASets

%% equation breaker
function [check,isterminal,direction] = EventsFcn(t,y,tstart)
% This function checks to see whether an ODE solver has exceeded a time
% limit to solve. The avoids numerical problems which can occur in systems
% with parameters and abundances on vastly different orders of magnitude.
global EventOccur
check = toc(tstart)<300; % The value that we want to be zero
if check == 0
    EventOccur = 1; % update that integration took too long
    disp('too long')
end
isterminal = 1;  % Halt integration
direction = 0;   % The zero can be approached from either direction
end

%% ode system
function eqns = odesys(t,u,k,MaxTot)
eqns = zeros(10,1); % To start with we have twelve empty equations
% Using u = [C  Cs  B Bs  I E1 C1 C2 C3 C4]
% Using k = [a1,d1,k1,a2,d2,k2,a3,d3,k3,a4,d4,k4]
u(u<-1) = 0; % check that abundances never go below zero
u(u>MaxTot) = MaxTot; % rough check that abundances do not exceed the largest maximum
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