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
global microwave_control;
global dt;
global probability;

% Create the system
system = spine.systems.system_1_spin();
system.setLarmorFrequency(2 * pi * 1e9);
system.setRabiFrequency(2 * pi * 10e6);

% Generate the driving signal
dt = 10e-12;
tpi = pi / system.getRabiFrequency();
Npi = round(tpi / dt);
t = (0:Npi-1) * dt;
microwave_control = cos(system.getLarmorFrequency() * (t - 0.5*dt));

% Simulate the system
nsim = 0;
Nsim = Npi;
spine.simulate(system.getDimension(), @inHamiltonian, @outOperation, @spine.solvers.solver_analytical_xz);

% Plot the result
figure();
plot(t, probability);

function [run, H, timestep] = inHamiltonian()
    global nsim;
    global Nsim;
    global system;
    global microwave_control;
    global dt;

    if (nsim < Nsim)
        
        % Set the signal at this time instance
        system.setMicrowaveControl(microwave_control(nsim + 1));
        
        % Update the Hamiltonian accordingly
        H = system.updateHamiltonian();
        
        % Provide the current timestep, and continue the simulation
        timestep = dt;
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
    global dt;
    global probability;

    % Plot every Nplot-th point
    Nplot = 5;
    
    % Store the measurement probability
    probability(nsim + 1) = abs(U(1,1))^2;
    
    % Show the rotation in the Bloch sphere
    if (mod(nsim, Nplot) == 0)
        system.plot(U, nsim * dt, 2);
        pause(0.001);
    end
    
    % Continue the simulation
    nsim = nsim + 1;
    if (nsim == Nsim)
        
        % Last point, print the fidelity and unitary matrix
        F = spine.fidelity(system.getDimension(), U, pi, 0);
        disp(F);
        disp(U);
        
    end
    
end
