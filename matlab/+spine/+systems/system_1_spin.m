%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef system_1_spin < spine.systems.system_spin

    methods
        
        function obj = system_1_spin()
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
            H = [-obj.larmor_frequency / 2.0, obj.rabi_frequency * obj.microwave_control;
                 obj.rabi_frequency * obj.microwave_control, obj.larmor_frequency / 2.0];
        end

    end
    
end
