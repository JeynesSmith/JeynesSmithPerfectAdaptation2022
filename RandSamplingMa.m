% This script is used to generate data regarding how parameters affect the
% ability for the Michaelian (Ma et al.) network to generate the RPA response,
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
maxI = 1e15; % max accepted input
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
    init = Parameters(i,9:10);
    kinit = ones(1,12);
    kinit(1:8) = Parameters(i,1:8);
    kinit(10:11) = Parameters(i,9:10);
    kinit(9) = Parameters(i,11);
    theoreticalEndPoint = kinit(4)/kinit(3)*kinit(9); % define the estimated setpoint
    
    % calculate first input for which RPA can be observed
    startI = 1e-5; % initialise input
    tempAbund = zeros(2,2); % initialise abundance matrix
    % broad search (increasing by order of magnitude at a time)
    for j = 1:2
        k = [kinit(1:8) startI/changeInput(j) kinit(9:end)];
        tstart=tic;
        [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, init, odeset('Events',@(t,y) EventsFcn(t,y,tstart))); % simulate system
        tempAbund(j,:) = u(end,:);
    end
    while ~(abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ... % check if A* within tolerance of setpoint
            && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ... % check if A* changed within tolerance
            && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2))>=Btol) ... % check if B* changed more than tolerance
            && maxI > startI % stop if input gets too large
        startI = startI*10; % increase input by order of magnitude
        tempAbund(2,:) = tempAbund(1,:);
        k = [kinit(1:8) startI kinit(9:end)];
        tstart=tic; % initialise timer for time check
        [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, init, odeset('Events',@(t,y) EventsFcn(t,y,tstart))); % simulate system
        tempAbund(1,:) = u(end,:); % update abundance
    end
    % refined search (increasing by order of 0.1*startI from broad search)
    startI = startI/5; % decrease input by one order of magnitude and then increase
    stepSize = startI/2; % set step size
    count = 2; % initialise count
    for j = 1:2
        k = [kinit(1:8) startI/j kinit(9:end)];
        tstart=tic;
        [t,u] = ode23s(@(t,u) odesys(t,u,k), tspan, init, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
        tempAbund(j,:) = u(end,:);
    end
    % repeat checks and updates for smaller input increase
    while ~(abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
            && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
            && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2))>=Btol) ...
            && maxI > startI && count<maxCount % added check to see if returned to previous order of magnitude
        startI = startI + stepSize; % increase by stepsize
        tempAbund(2,:) = tempAbund(1,:);
        k = [kinit(1:8) startI kinit(9:end)];
        tstart=tic;
        [t,u] = ode23s(@(t,u) odesys(t,u,k), tspan, init, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
        tempAbund(1,:) = u(end,:);
        count = count + 1;
    end
    
    % if RPA, calculate end input
    if maxI ~= startI % check that a starting input was found
        
        finalInput = startI;  % set the final input to be the start input
        % repeat updates until we reach an order of magnitude or until RPA is lost
        while abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
                && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
                && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2))>Btol ...
                && maxI > finalInput && count<maxCount
            finalInput = finalInput + stepSize;
            tempAbund(2,:) = tempAbund(1,:);
            k = [kinit(1:8) finalInput kinit(9:end)];
            tstart=tic;
            [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, init, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
            tempAbund(1,:) = u(end,:);
            count = count + 1;
        end
        
        % check that we still have RPA
        if abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint)<=Atol ...
                && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
                && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2))>Btol
            
            tempAbund(2,:) = tempAbund(1,:);
            finalInput = finalInput*10;
            k = [kinit(1:8) finalInput kinit(9:end)];
            tstart=tic;
            [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, init, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
            tempAbund(1,:) = u(end,:);
            simNum = 1;
            % broad search (increasing by order of magnitude at a time)
            % repeat checks and updates
            while abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint) <= Atol ...
                    && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
                    && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2)) > Btol ...
                    && maxSim >= simNum
                finalInput = finalInput*10;
                tempAbund(2,:) = tempAbund(1,:);
                k = [kinit(1:8) finalInput kinit(9:end)];
                tstart=tic;
                [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, init, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
                tempAbund(1,:) = u(end,:);
                simNum = simNum + 1;
            end
            % refined search (increasing by order of 0.1*startI from broad search)
            stepSize = finalInput/10; % set step size
            finalInput = max(startI,finalInput/5); % decrease input by one order of magnitude and then increase
            count = 2;
            for j = 1:2
                k = [kinit(1:8) finalInput/j kinit(9:end)];
                tstart=tic;
                [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, init, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
                tempAbund(j,:) = u(end,:);
            end
            % repeat checks and updates for smaller stepsize
            while abs(tempAbund(1,1) - theoreticalEndPoint)/abs(theoreticalEndPoint) <= Atol ...
                    && abs(tempAbund(1,1) - tempAbund(2,1))/abs(tempAbund(1,1))<=2*Atol ...
                    && abs(tempAbund(1,2) - tempAbund(2,2))/abs(tempAbund(1,2)) > Btol ...
                    && count<maxCount
                finalInput = finalInput + stepSize;
                tempAbund(2,:) = tempAbund(1,:);
                k = [kinit(1:8) finalInput kinit(9:end)];
                tstart=tic;
                [~,u] = ode23s(@(t,u) odesys(t,u,k), tspan, init, odeset('Events',@(t,y) EventsFcn(t,y,tstart)));
                tempAbund(1,:) = u(end,:);
                simNum = simNum + 1;
                count = count + 1;
            end
        end
        if EventOccur == 1 % check if simulation ever took too long to run
            disp('had to break')
            RPASets(i,:) = [Parameters(i,:) -1 -1];  % save in RPASets but with negative values of Is and If
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
%% save the parameter set 
save MaRSTot nonRPASets RPASets

%% equation breaker
function [position,isterminal,direction] = EventsFcn(t,y,tstart)
global EventOccur
position = toc(tstart)<300; % The value that we want to be zero
if position == 0
    EventOccur = 1;
    disp('too long')
end
isterminal = 1;  % Halt integration
direction = 0;   % The zero can be approached from either direction
end

%% ode system
function eqns = odesys(t,u,k)
eqns = zeros(2,1); % To start with we have twelve empty equations
% Using u = [A B]
% Using k = [k1,k2,k3,k4,K1,K2,K3,K4,I,E1,At,Bt]
u(u<0) = 0; % check that abundances never go below zero
temp = k(11:12)'; % get total substrate abundances
u(u>temp) = temp(u>temp); % check that abundances never go above total substrate abundances
eqns(1) = k(1)*k(9)*(k(11)-u(1))/(k(5) + k(11)-u(1)) - k(2)*u(1)*u(2)/(k(6) + u(1));
eqns(2) = k(3)*u(1)*(k(12)-u(2))/(k(7) + k(12)-u(2)) - k(4)*k(10)*u(2)/(k(8) + u(2));
end