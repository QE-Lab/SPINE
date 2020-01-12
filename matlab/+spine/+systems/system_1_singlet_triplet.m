%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef system_1_singlet_triplet < handle
    
    properties
        
        % Hamiltonian Properties
        magnetic_gradient = 0;
        
        % Hamiltonian Control
        exchange_interaction = 0;
        
        % Plot handles
        p = gobjects(1, 1);
        
    end
    
    methods
        
        function obj = system_1_singlet_triplet()
            
        end
        
        function dimension = getDimension(~)
            dimension = 2;
        end
        
        function state = initialize(~)
            state = [1; 0];
        end
        
        function [Px, Py, Pz] = measure(~, state)
            Px = abs((state(1) + 1 * state(2))/sqrt(2)).^2;
            Py = abs((state(1) + 1i * state(2))/sqrt(2)).^2;
            Pz = abs(state(1)).^2;
        end
        
        function H = updateHamiltonian(obj)
            H = [obj.exchange_interaction / 2.0, obj.magnetic_gradient / 2.0;
                 obj.magnetic_gradient / 2.0, -obj.exchange_interaction / 2.0];
        end
        
        % Hamiltonian Property Getters/Setters
        function setMagneticGradient(obj, magnetic_gradient)
            obj.magnetic_gradient = magnetic_gradient;
        end
        
        function magnetic_gradient = getMagneticGradient(obj)
            magnetic_gradient = obj.magnetic_gradient;
        end
        
        % Hamiltonian Control Setters
        function setExchangeInteraction(obj, exchange_interaction)
            obj.exchange_interaction = exchange_interaction;
        end
        
        % Plot
        function plot(obj, state_or_U)

            % Init, rotate, measure
            if (size(state_or_U, 2) == 1)
                state = state_or_U;
            else
                state = state_or_U * obj.initialize();
            end
            [Px, Py, Pz] = obj.measure(state);

            % Plot
            if (isempty(fieldnames(get(obj.p))))
                figure();
                spine.plotBlochSphere();
                obj.p(1) = plot3(2*Px-1, 2*Py-1, 2*Pz-1);
            else
                obj.p(1).XData = [obj.p(1).XData, 2*Px-1];
                obj.p(1).YData = [obj.p(1).YData, 2*Py-1];
                obj.p(1).ZData = [obj.p(1).ZData, 2*Pz-1];
            end
    
        end
        
    end
    
end
