%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef system_dispersive_readout < handle
    
    properties
        
        % Boltzmann constant (eV)
        kB = 8.62e-5;
        
        % Hamiltonian Properties
        electron_temperature = 0;
        tunnel_rate = 0;
        lever_arm = 0;
        
        % Hamiltonian Control
        gate_voltage = 0;
        
        % Plot handles
        p = gobjects(1, 1);
        
    end
    
    methods
        
        function obj = system_dispersive_readout()
            
        end
        
        function dimension = getDimension(obj)
            dimension = 2;
        end
        
        function state = initialize(~, P)
            state = [1; P];
        end
        
        function [P] = measure(~, state)
            P = state(2);
        end
        
        function H = updateHamiltonian(obj)
            H = 1i * [0, 0; obj.tunnel_rate / (1 + exp(obj.lever_arm * obj.gate_voltage/obj.kB/obj.electron_temperature)), -obj.tunnel_rate];
        end
        
        % Hamiltonian Property Getters/Setters
        function setElectronTemperature(obj, electron_temperature)
            obj.electron_temperature = electron_temperature;
        end
        
        function electron_temperature = getElectronTemperature(obj)
            electron_temperature = obj.electron_temperature;
        end
        
        function setTunnelRate(obj, tunnel_rate)
            obj.tunnel_rate = tunnel_rate;
        end
        
        function tunnel_rate = getTunnelRate(obj)
            tunnel_rate = obj.tunnel_rate;
        end
        
        function setLeverArm(obj, lever_arm)
            obj.lever_arm = lever_arm;
        end
        
        function lever_arm = getLeverArm(obj)
            lever_arm = obj.lever_arm;
        end
        
        % Hamiltonian Control Setters
        function setGateVoltage(obj, gate_voltage)
            obj.gate_voltage = gate_voltage;
        end
        
        % Plot
        function plot(obj, state_or_U, t)

            % Init, rotate, measure
             if (size(state_or_U, 2) == 1)
                state = state_or_U;
            else
                state = state_or_U * obj.initialize(0);
             end
            [P] = obj.measure(state);

            % Plot
            if (isempty(fieldnames(get(obj.p))))
                figure();
                obj.p(1) = plot(t, P);
            else
                obj.p(1).XData = [obj.p(1).XData, t];
                obj.p(1).YData = [obj.p(1).YData, P];
            end
    
        end
        
    end
    
end
