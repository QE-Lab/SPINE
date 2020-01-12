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
global microwave_amplitude;
global microwave_phase;
global dt;
global probability;

% Create the system
system = spine.systems.system_1_spin_rwa();
system.setMicrowaveFrequency(2 * pi * 1e9);
system.setLarmorFrequency(2 * pi * 1e9);
system.setRabiFrequency(2 * pi * 1e6);

% Generate the driving signal
dt = 50e-12;
tpi = pi / system.getRabiFrequency();
Npi = tpi / dt;
t = (1:Npi) * dt;
microwave_amplitude = ones(1, Npi);
microwave_phase = [zeros(1, Npi/2), pi/2*ones(1, Npi/2)];

% Simulate the system
nsim = 0;
Nsim = Npi;
spine.simulate(system.getDimension(), @inHamiltonian, @outOperation, @spine.solvers.solver_expm);

% Plot the result
figure();
plot(t, probability);

function [run, H, timestep] = inHamiltonian()
    global nsim;
    global Nsim;
    global system;
    global microwave_amplitude;
    global microwave_phase;
    global dt;

    if (nsim < Nsim)
        
        % Set the signal at this time instance
        system.setMicrowaveAmplitude(microwave_amplitude(nsim + 1));
        system.setMicrowavePhase(microwave_phase(nsim + 1));
        
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
    global probability;
    
    % Plot every Nplot-th point
    Nplot = 5;
    
    % Store the measurement probability
    probability(nsim + 1) = abs(U(1,1))^2;
    
    % Show the rotation in the Bloch sphere
    if (mod(nsim, Nplot) == 0)
        system.plot(U, 0, 1);  % t = 0, as we are already in the rotating frame!
        pause(0.001);
    end
    
    % Continue the simulation
    nsim = nsim + 1;
    if (nsim == Nsim)
        
        % The ideal X, pi/2 rotation
        Ux = [cos(pi/4), -1i * (cos(0) + 1i * sin(0)) * sin(pi/4);
              -1i * (cos(0) - 1i * sin(0)) * sin(pi/4), cos(pi/4)];
        Uy = [cos(pi/4), -1i * (cos(pi/2) + 1i * sin(pi/2)) * sin(pi/4);
              -1i * (cos(pi/2) - 1i * sin(pi/2)) * sin(pi/4), cos(pi/4)];
        Uideal = Uy * Ux;
        
        % Last point, print the fidelity and unitary matrix
        F = spine.fidelity(system.getDimension(), U, Uideal);
        disp(F);
        disp(U);
        
    end
    
end
