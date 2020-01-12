%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef system_1_spin_rwa < spine.systems.system_spin
    
    properties

        % Hamiltonian Properties
        microwave_frequency = 0;
        
        % Hamiltonian Control
        microwave_amplitude = 0;
        microwave_phase = 0;

    end
    
    methods
        
        function obj = system_1_spin_rwa()
            obj = obj@spine.systems.system_spin(1);
        end
        
        function index = getIndex(~, states)
            if (states(1) == spine.systems.system_spin.STATE_0)
                index = 1;
            elseif (states(1) == spine.systems.system_spin.STATE_1)
                index = 2;
            end
        end
        
        function index = getIndexMeasurement(~, dot, state)
            index = [];
            if (dot == 1)
                if (state == spine.systems.system_spin.STATE_0)
                    index = 1;
                elseif (state == spine.systems.system_spin.STATE_1)
                    index = 2;
                end
            end
        end
        
        function dimension = getDimension(obj)
            dimension = 2;
        end

        function H = updateHamiltonian(obj)
            H = [(obj.microwave_frequency - obj.larmor_frequency) / 2.0, obj.rabi_frequency * obj.microwave_amplitude / 2.0 * (cos(obj.microwave_phase) + 1i * sin(obj.microwave_phase));
                 obj.rabi_frequency * obj.microwave_amplitude / 2.0 * (cos(obj.microwave_phase) - 1i * sin(obj.microwave_phase)), -(obj.microwave_frequency - obj.larmor_frequency) / 2.0];
        end
        
        % Hamiltonian Property Getters/Setters
        function setMicrowaveFrequency(obj, value)
            obj.microwave_frequency = value;
        end

        function microwave_frequency = getMicrowaveFrequency(obj)
            microwave_frequency = obj.microwave_frequency;
        end

        % Hamiltonian Control Setters
        function setMicrowaveControl(~, varargin)
            error('Use setMicrowaveAmplitude() and setMicrowavePhase() instead!');
        end
        
        function setMicrowaveAmplitude(obj, microwave_amplitude)
            obj.microwave_amplitude = microwave_amplitude;
        end
        
        function setMicrowavePhase(obj, microwave_phase)
            obj.microwave_phase = microwave_phase;
        end
        
    end
    
end
