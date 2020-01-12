%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef system_2_spin_1_singlet < spine.systems.system_spin

    methods
        
        function obj = system_2_spin_1_singlet()
            obj = obj@spine.systems.system_spin(2);
        end
        
        function index = getIndex(~, states)
            if ((states(1) == spine.systems.system_spin.STATE_0) && (states(2) == spine.systems.system_spin.STATE_0))
                index = 1;
            elseif ((states(1) == spine.systems.system_spin.STATE_0) && (states(2) == spine.systems.system_spin.STATE_1))
                index = 2;
            elseif ((states(1) == spine.systems.system_spin.STATE_1) && (states(2) == spine.systems.system_spin.STATE_0))
                index = 3;
            elseif ((states(1) == spine.systems.system_spin.STATE_1) && (states(2) == spine.systems.system_spin.STATE_1))
                index = 4;
            elseif ((states(1) == spine.systems.system_spin.STATE_N) && (states(2) == spine.systems.system_spin.STATE_S))
                index = 5;
            end
        end
        
        function index = getIndexMeasurement(~, dot, state)
            index = [];
            if (dot == 1)
                if (state == spine.systems.system_spin.STATE_0)
                    index = [1, 2];
                elseif (state == spine.systems.system_spin.STATE_1)
                    index = [3, 4];
                elseif (state == spine.systems.system_spin.STATE_N)
                    index = 5;
                end
            elseif (dot == 2)
                if (state == spine.systems.system_spin.STATE_0)
                    index = [1, 3];
                elseif (state == spine.systems.system_spin.STATE_1)
                    index = [2, 4];
                elseif (state == spine.systems.system_spin.STATE_S)
                    index = 5;
                end
            end
        end

        function dimension = getDimension(obj)
            dimension = 5;
        end
        
        function H = updateHamiltonian(obj)
            
            % Single microwave driveline
            microwave_control = obj.microwave_control(1) + obj.microwave_control(2);
            
            % Common Rabi frequency in simple Hamiltonian
            rabi_frequency = obj.rabi_frequency(1);
            
            % All tunnel rates are the same for a double dot
            tunnel_control = obj.tunnel_control(1, 1);
            
            % Only singlet for dot 2
            charging_energy = obj.charging_energy(2);
            detuning_control = obj.detuning_control(2);
            
            H = [-(obj.larmor_frequency(1) + obj.larmor_frequency(2)) / 2.0, microwave_control*rabi_frequency, microwave_control*rabi_frequency, 0, 0;
                microwave_control*rabi_frequency, -(obj.larmor_frequency(1) - obj.larmor_frequency(2)) / 2.0, 0, microwave_control*rabi_frequency, tunnel_control;
                microwave_control*rabi_frequency, 0, +(obj.larmor_frequency(1) - obj.larmor_frequency(2)) / 2.0, microwave_control*rabi_frequency, -tunnel_control;
                0, microwave_control*rabi_frequency, microwave_control*rabi_frequency, +(obj.larmor_frequency(1) + obj.larmor_frequency(2)) / 2.0, 0;
                0, tunnel_control, -tunnel_control, 0, charging_energy - detuning_control];
           
        end
        
    end
    
end