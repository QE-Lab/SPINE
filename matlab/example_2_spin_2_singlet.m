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
global dt;
global probability;

% Create the system
system = spine.systems.system_2_spin_2_singlet();
system.setChargingEnergy(2 * pi * 100e9);
system.setRabiFrequency(2 * pi * 1e6);
system.setTunnelControl(2 * pi * 10e6 * sqrt(2));
system.setLarmorFrequency(1, 2 * pi * 1e9);
system.setLarmorFrequency(2, 2 * pi * 1e9);

% Generate the driving signal
dt = 50e-12;
Nop = round(500e-9 / dt);
t = (1:Nop) * dt;
system.setDetuningControl(1, (1 - 2e-3) * system.getChargingEnergy());

% Plot the eigen energies
figure();
spine.plotEigenEnergies(system);

% Simulate the system
nsim = 0;
Nsim = Nop;
state = [1; 1; 0; 0; 0; 0] / sqrt(2);
spine.simulate(system.getDimension(), @inHamiltonian, @outOperation, @spine.solvers.solver_expm, state);

% Plot the result
figure();
plot(t, probability);

function [run, H, timestep] = inHamiltonian()
    global nsim;
    global Nsim;
    global system;
    global dt;

    if (nsim < Nsim)
        
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
    Nplot = 10;
    
    % Store the measurement probability of the singlet states
    probability(nsim + 1) = abs(U(5,1))^2 + abs(U(6,1))^2;
    
    % Show the rotation in the Bloch sphere
    if (mod(nsim, Nplot) == 0)
        system.plot(U, nsim * dt);
        pause(0.001);
    end
    
    % Continue the simulation
    nsim = nsim + 1;
    if (nsim == Nsim)
        
        % Last point, print the unitary matrix
        disp(U);
        
    end
    
end
