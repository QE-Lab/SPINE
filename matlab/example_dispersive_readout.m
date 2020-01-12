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
global gate_voltage;
global dt;
global probability;

% Create the system
system = spine.systems.system_dispersive_readout();
system.setElectronTemperature(0.1);
system.setTunnelRate(2 * pi * 1e9);
system.setLeverArm(0.05);

% Generate the driving signal
% Drive frequency determines waveform shape (capacitive vs resistive regimes)
w0 = 2 * system.getTunnelRate();
tmax = 10 * pi / w0;
Npts = 100;
dt = tmax / Npts;
t = (1:Npts) * dt;
gate_voltage = 0.4e-3 / system.getLeverArm() * sin(w0 * t);

% Simulate the system
nsim = 0;
Nsim = Npts;
spine.simulate(system.getDimension(), @inHamiltonian, @outOperation, @spine.solvers.solver_expm);

% Plot the result
figure();
subplot(2, 1, 1);
plot(t(2:end), gate_voltage(2:end));
ylabel('Excitation voltage');
subplot(2, 1, 2);
q = 1.602e-19;
I = -q * diff(probability) ./ diff(t);
plot(t(2:end), I);
ylabel('Gate current');

function [run, H, timestep] = inHamiltonian()
    global nsim;
    global Nsim;
    global system;
    global gate_voltage;
    global dt;

    if (nsim < Nsim)
        
        % Set the signal at this time instance
        system.setGateVoltage(gate_voltage(nsim + 1));
        
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
    global system;
    global dt;
    global probability;
    
    % Store the measurement probability
    probability(nsim + 1) = U(2,1);
    
    % Plot every Nplot-th point
    Nplot = 1;
    
    % Show the rotation in the Bloch sphere
    if (mod(nsim, Nplot) == 0)
        system.plot(U, nsim * dt);
        pause(0.001);
    end
    
    % Continue the simulation
    nsim = nsim + 1;
    
end
