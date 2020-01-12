%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cleanup
clear all;
close all;
clc;

% Global variables
global nsim;
global Nsim;
global system;
global J t;
global probability;

% Create the system
system = spine.systems.system_1_singlet_triplet();
system.setMagneticGradient(160e6);

% Generate the driving signal, consisting of two levels with different duration
J1 = system.getMagneticGradient() * (1 + sqrt(2));
J2 = system.getMagneticGradient() * (J1 - system.getMagneticGradient()) / (J1 + system.getMagneticGradient());
J = [J1, J2];
t = pi ./ sqrt(system.getMagneticGradient()^2 + J.^2);

% Simulate the system
nsim = 0;
Nsim = 50;
spine.simulate(system.getDimension(), @inHamiltonian, @outOperation, @spine.solvers.solver_analytical_xz);

% Plot the result
figure();
t1 = (1:Nsim) * t(1) / Nsim;
t2 = (1:Nsim) * t(2) / Nsim;
plot([t1, t(1) + t2], probability, '.');

function [run, H, timestep] = inHamiltonian()
    global nsim;
    global Nsim;
    global system;
    global J t;

    if (nsim < Nsim)
        
        % Set the signal at this time instance
        system.setExchangeInteraction(J(1));
        
        % Update the Hamiltonian accordingly
        H = system.updateHamiltonian();
        
        % Provide the current timestep, and continue the simulation
        timestep = t(1) / Nsim;
        run = true;
        
    elseif (nsim < 2*Nsim)
        
        % Set the signal at this time instance
        system.setExchangeInteraction(J(2));
        
        % Update the Hamiltonian accordingly
        H = system.updateHamiltonian();
        
        % Provide the current timestep, and continue the simulation
        timestep = t(2) / Nsim;
        run = true;
        
    else
        run = false;
        H = [];
        timestep = [];
    end
    
end

function outOperation(U)
    global nsim;
    global Nsim;
    global system;
    global probability;
    
    % Plot every Nplot-th point
    Nplot = 5;

    % Store the measurement probability
    probability(nsim + 1) = abs(U(1,1))^2;
    
    % Show the rotation in the Bloch sphere
    if (mod(nsim, Nplot) == 0)
        system.plot(U);
        pause(0.001);
    end
    
    % Continue the simulation
    nsim = nsim + 1;
    if (nsim == 2*Nsim)
        
        % Last point, print the fidelity and unitary matrix
        F = spine.fidelity(system.getDimension(), U, pi/2, -pi/2);
        disp(F);
        disp(U);
        
    end
    
end
