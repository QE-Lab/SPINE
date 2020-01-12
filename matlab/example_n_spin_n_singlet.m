%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cleanup
clear all;
close all;
clc;

% Global variables
global nsim;
global system;
global dt;
global N1 N2 N3 N4 N5 N6 N7 N8 N9;
global clk_I1 clk_I3 clk_Q3 clk_I4 clk_Q4;
global tunnel_coupling_swap tunnel_coupling_cz tunnel_coupling_measure

% Create the system
n_dots = 4;
system = spine.systems.system_n_spin_n_singlet(n_dots);
for i=1:n_dots
    system.setChargingEnergy(i, 2 * pi * 100e9);
    system.setDetuningControl(i, 2 * pi * 0e9)
    system.setLarmorFrequency(i, 2 * pi * 2.5e9);
    system.setRabiFrequency(i, 2 * pi * 10e6);
end
system.setDetuningControl(3, 2 * pi * 75e9);
system.setDetuningControl(4, 2 * pi * 75e9);
system.setLarmorFrequency(3, 2 * pi * 2.5e9+0.5e9*sqrt(2));
system.setLarmorFrequency(4, 2 * pi * 2.5e9-0.5e9*sqrt(2));

tunnel_coupling_swap = 2*pi*0.5e9 * 1.6;
tunnel_coupling_cz = 2*pi*0.5e9 * 1.6;
tunnel_coupling_measure = 2*pi*0.5e9 * 100;
    
% The system is tuned such that each operation has a speed of 10 MHz
% omega_SWAP = 2*pi*10e6;
% omega_CZ = 2*pi*10e6;
omega_R = 2*pi*10e6;
    
% Simulation time step
dt = 50e-12;

% Perform the following simulation
Npi = round(pi/omega_R/dt);
N0 = 0;             % -            initialize the system to the ground state
N1 = N0 + Npi;      % 50.0 ns      rotate q1 from 0 to 1 using a microwave signal
N2 = N1 + Npi;      % 50.0 ns      swap q1 and q2 using the tunnel coupling
N3 = N2 + Npi/2;    % 25.0 ns      rotate q3 and q4 in plane using a microwave signal
N4 = N3 + Npi/4;    % 12.5 ns      adiabatic change the tunnel coupling between q1/q3 and q2/q4
N5 = N4 + Npi/2;    % 25.0 ns      controlled-z between q1/q3 and q2/q4
N6 = N5 + Npi/4;    % 12.5 ns      adiabatic change the tunnel coupling between q1/q3 and q2/q4
N7 = N6 + Npi/2;    % 25.0 ns      rotate q3 and q4 out plane using a microwave signal, out of phase
N8 = N7 + Npi;      % 50.0 ns      measure q3/q4
N9 = N8 + Npi;      % 50.0 ns      rotate q1 using a microwave signal
N = N9;
    
% Generate the I/Q clocks for the simulation
t = (1:N) * dt;
clk_I1 = sin(system.getLarmorFrequency(1) * t);
clk_I3 = sin(system.getLarmorFrequency(3) * t);
clk_Q3 = cos(system.getLarmorFrequency(3) * t);
clk_I4 = sin(system.getLarmorFrequency(4) * t);
clk_Q4 = cos(system.getLarmorFrequency(4) * t);

% -------------------------- INITIALIZE -------------------------------
state = system.initialize();

% Simulate the system
nsim = 0;
tic;
spine.simulate(system.getDimension(), @inHamiltonian, @outOperation, @spine.solvers.solver_taylor_sparse_approx, state);
toc;

function [run, H, timestep] = inHamiltonian()
    global nsim;
    global N1 N2 N3 N4 N5 N6 N7 N8 N9;
    global clk_I1 clk_I3 clk_Q3 clk_I4 clk_Q4;
    global tunnel_coupling_swap tunnel_coupling_cz tunnel_coupling_measure
    global system;
    global dt;
    global t0;

    if (nsim < N1)
        
        % Apply the microwave signal to q1
        system.setMicrowaveControl(1, clk_I1(nsim+1));
        
    elseif (nsim < N2)
        
        % Change the tunnel coupling diabatically
        system.setTunnelControl(1, 2, tunnel_coupling_swap);
        
    elseif (nsim < N3)
        
        % Apply the microwave signal to q3 and q4
        system.setTunnelControl(1, 2, 0);
        system.setMicrowaveControl(3, clk_I3(nsim+1));
        system.setMicrowaveControl(4, clk_I4(nsim+1));
        t0 = 0;
        
    elseif (nsim < N4)
        
        % Ramp the tunnel coupling
        dt0 = tunnel_coupling_cz / (N4 - N3);
        t0 = t0 + dt0;
        system.setTunnelControl(1, 3, t0);
        system.setTunnelControl(2, 4, t0);
        
    elseif (nsim < N5)
        
        % Hold the tunnel coupling
        system.setTunnelControl(1, 3, t0);
        system.setTunnelControl(2, 4, t0);
        
    elseif (nsim < N6)
        
        % Ramp the tunnel coupling
        dt0 = tunnel_coupling_cz / (N4 - N3);
        t0 = t0 - dt0;
        system.setTunnelControl(1, 3, t0);
        system.setTunnelControl(2, 4, t0);
        
    elseif (nsim < N7)
        
        % Apply the microwave signal to q3 and q4
        system.setMicrowaveControl(3, clk_Q3(nsim+1));
        system.setMicrowaveControl(4, clk_Q4(nsim+1));
        t0 = 0;
        
    elseif (nsim < N8)
        
        % Detune and ramp the tunnel coupling
        system.setDetuningControl(3, 2*pi*0e9);
        system.setDetuningControl(4, 2*pi*105e9);
        dt0 = tunnel_coupling_measure / (N8 - N7);
        t0 = t0 + dt0;
        system.setTunnelControl(3, 4, t0);
        
    elseif (nsim < N9)
        
        % Apply the microwave signal to q1
        system.setMicrowaveControl(1, -clk_I1(nsim+1));

    else
        run = false;
        H = [];
        timestep = [];
        return
    end
    
    % Update the Hamiltonian accordingly
    H = system.updateHamiltonian();

    % Provide the current timestep, and continue the simulation
    timestep = dt;
    run = true;
    
end

function outOperation(U)
    global nsim;
    global dt;
    global system;
    global N9;

    % Plot every Nplot-th point
    Nplot = 10;
    
    % Show the rotation in the Bloch sphere
    if (mod(nsim, Nplot) == 0)
        system.plot(U, nsim * dt, 1, 2);
        pause(0.001);
    end
    
    % Continue the simulation
    nsim = nsim + 1;
    if (nsim == N9)
        
        % Last point, print the state
        disp(U);
        
    end
    
end

