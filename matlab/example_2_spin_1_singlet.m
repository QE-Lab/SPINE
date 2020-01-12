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
global microwave_control_1 microwave_control_2 detuning_control;
global dt;
global probability;

% Create the system
system = spine.systems.system_2_spin_1_singlet();
system.setChargingEnergy(2 * pi * 100e9);
system.setRabiFrequency(2 * pi * 1e6);
system.setTunnelControl(2 * pi * 10e6 * sqrt(2));
system.setLarmorFrequency(1, 2 * pi * (1e9 - 10e6));
system.setLarmorFrequency(2, 2 * pi * (1e9 + 10e6));

% Generate the driving signal
dt = 50e-12;
Ngauss = round(pi/2/system.getRabiFrequency()/dt);
ygauss = exp(-((1:Ngauss) - 0.5*Ngauss).^2/2/(0.15*Ngauss).^2);
ygauss = ygauss * Ngauss/sum(ygauss);
yzeros = zeros(1, Ngauss);
ydetun = ones(1, Ngauss) * (1 - 2e-3) * system.getChargingEnergy();
ypiby2 = ones(1, Ngauss) * pi/2;

% The microwave signals
Nop = 7 * Ngauss;
t = (1:Nop) * dt;
microwave_control_1 = [ygauss, yzeros, ygauss, yzeros, ygauss, yzeros, ygauss] .* cos(system.getLarmorFrequency(1) * t ...
                    + [yzeros, yzeros, ypiby2, yzeros, yzeros, yzeros, ypiby2]);
microwave_control_2 = [yzeros, yzeros, yzeros, 2*ygauss, yzeros, yzeros, yzeros] .* cos(system.getLarmorFrequency(2) * t ...
                    + [yzeros, yzeros, yzeros, yzeros, yzeros, yzeros, yzeros]);
detuning_control    = [yzeros, ydetun, yzeros, yzeros, yzeros, ydetun, yzeros];

% Plot the eigen energies
figure();
spine.plotEigenEnergies(system);

% Simulate the system
nsim = 0;
Nsim = Nop;
spine.simulate(system.getDimension(), @inHamiltonian, @outOperation, @spine.solvers.solver_expm);

% Plot the result
figure();
plot(t, probability);

function [run, H, timestep] = inHamiltonian()
    global nsim;
    global Nsim;
    global system;
    global microwave_control_1 microwave_control_2 detuning_control;
    global dt;

    if (nsim < Nsim)
        
        % Set the signal at this time instance
        system.setMicrowaveControl(1, microwave_control_1(nsim + 1));
        system.setMicrowaveControl(2, microwave_control_2(nsim + 1));
        system.setDetuningControl(detuning_control(nsim + 1));
        
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
    
    % Store the measurement probability of the singlet state
    probability(nsim + 1) = abs(U(5,1))^2;
    
    % Show the rotation in the Bloch sphere
    if (mod(nsim, Nplot) == 0)
        system.plot(U, nsim * dt, 1);
        pause(0.001);
    end
    
    % Continue the simulation
    nsim = nsim + 1;
    if (nsim == Nsim)
        
        % Last point, print the unitary matrix
        disp(U);
        
    end
    
end
