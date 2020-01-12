%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef system_spin < handle

        properties
            
            % Number of dots
            dots;
            
            % Hamiltonian Properties
            larmor_frequency;
            rabi_frequency;
            charging_energy;
            singlet_triplet_energy;

            % Hamiltonian Control
            microwave_control;
            detuning_control;
            tunnel_control;
            
            % Plot handles
            p;
            
        end
        
        properties (Constant)
            
            % Supported spin states:
            STATE_N = 0              % 0 electrons
            STATE_0 = 1              % 1 electron in ground state
            STATE_1 = 2              % 1 electron in excited state
            STATE_S = 3              % 2 electrons in singlet configuration (or in 00 spin state, depending on the Hamiltonian)
            STATE_T0 = 4             % 2 electrons in triplet (0) configuration (or in 01 spin state, depending on the Hamiltonian)
            STATE_TP = 5             % 2 electrons in triplet (+) configuration (or in 10 spin state, depending on the Hamiltonian)
            STATE_TM = 6             % 2 electrons in triplet (-) configuration (or in 11 spin state, depending on the Hamiltonian)
            
        end
        
        methods (Abstract)
            getIndex(obj, states)
            getIndexMeasurement(obj, dot, state)
            getDimension(obj)
        end
        
        methods
           
            function obj = system_spin(dots)
                obj.dots = dots;
                obj.larmor_frequency = zeros(1, dots);
                obj.rabi_frequency = zeros(1, dots);
                obj.charging_energy = zeros(1, dots);
                obj.singlet_triplet_energy = zeros(1, dots);
                obj.microwave_control = zeros(1, dots);
                obj.detuning_control = zeros(1, dots);
                obj.tunnel_control = zeros(dots, dots);
                obj.p = gobjects(1, 2 * dots);
            end
            
            function state = initialize(obj)
                state = zeros(obj.getDimension(), 1);
                states = obj.STATE_0 * ones(1, obj.dots);
                state(obj.getIndex(states)) = 1;
            end
            
            function [Px, Py, Pz] = measure(obj, state, varargin)
                
                % Measure in rotating frame at time t if requested
                if (nargin == 3)
                    t = varargin{1};
                end

                % For every dot
                Px = zeros(1, obj.dots);
                Py = zeros(1, obj.dots);
                Pz = zeros(1, obj.dots);
                for dot=1:obj.dots
                    if (nargin == 3)
                        rot_x = cos(obj.larmor_frequency(dot) * t) + 1i * sin(obj.larmor_frequency(dot) * t);
                        rot_y = cos(obj.larmor_frequency(dot) * t + pi/2) + 1i * sin(obj.larmor_frequency(dot) * t + pi/2);
                    else
                        rot_x = 1;
                        rot_y = 1i;
                    end
                    
                    % Sum over all indices where this dot in the requested state
                    Px(dot) = sum(abs((state(obj.getIndexMeasurement(dot, obj.STATE_0)) + rot_x * state(obj.getIndexMeasurement(dot, obj.STATE_1)))/sqrt(2)).^2);
                    Py(dot) = sum(abs((state(obj.getIndexMeasurement(dot, obj.STATE_0)) + rot_y * state(obj.getIndexMeasurement(dot, obj.STATE_1)))/sqrt(2)).^2);
                    Pz(dot) = sum(abs(state(obj.getIndexMeasurement(dot, obj.STATE_0))).^2);
                end
            end
            
            function [Pn, Ps, Pt] = measureST(obj, state)
                
                % For every dot
                Pn = zeros(1, obj.dots);
                Ps = zeros(1, obj.dots);
                Pt = zeros(1, obj.dots);
                for dot=1:obj.dots
                    
                    % Sum over all indices where this dot in the requested state
                    Pn(dot) = sum(abs(state(obj.getIndexMeasurement(dot, obj.STATE_N))).^2);
                    Ps(dot) = sum(abs(state(obj.getIndexMeasurement(dot, obj.STATE_S))).^2);
                    Pt(dot) = sum(abs(state(obj.getIndexMeasurement(dot, obj.STATE_T0))).^2) + ...
                              sum(abs(state(obj.getIndexMeasurement(dot, obj.STATE_TP))).^2) + ...
                              sum(abs(state(obj.getIndexMeasurement(dot, obj.STATE_TM))).^2);
                end
            end
            
            function plot(obj, state_or_U, t, varargin)
                
                % Plot style lab
                % 0: No lab frame
                % 1: Arrow lab frame
                % 2: Trace lab frame
                plot_style_lab = 0;
                if (nargin == 4)
                    plot_style_lab = varargin{1};
                end
                
                % Plot singlet-triplet
                % 0: No additional plot
                % 1: S occupancy
                % 2: 1+S+T-N occupancy (expected number of electrons)
                plot_st = 0;
                if (nargin == 5)
                    plot_st = varargin{2};
                end
                
                % Init, rotate, measure
                if (size(state_or_U, 2) == 1)
                    state = state_or_U;
                else
                    state = state_or_U * obj.initialize();
                end
                if (plot_style_lab ~= 0)
                    [Px, Py, Pz] = obj.measure(state);
                end
                [Pxr, Pyr, Pzr] = obj.measure(state, t);
                if (plot_st ~= 0)
                    [Pn, Ps, Pt] = obj.measureST(state);
                end
                
                % Plot
                if (isempty(fieldnames(get(obj.p))))
                    figure();
                    for dot=1:obj.dots
                        subplot(floor(sqrt(obj.dots)), ceil(sqrt(obj.dots)), dot);
                        spine.plotBlochSphere();
                        if (plot_style_lab == 1)
                            obj.p(2*(dot-1)+1) = plot3([0, 2*Px(dot)-1], [0, 2*Py(dot)-1], [0, 2*Pz(dot)-1]);
                        elseif (plot_style_lab == 2)
                            obj.p(2*(dot-1)+1) = plot3(2*Px(dot)-1, 2*Py(dot)-1, 2*Pz(dot)-1);
                        else
                            obj.p(2*(dot-1)+1) = plot3(0, 0, 0);
                        end
                        obj.p(2*(dot-1)+2) = plot3(2*Pxr(dot)-1, 2*Pyr(dot)-1, 2*Pzr(dot)-1);
                    end
                    if (plot_st ~= 0)
                        figure();
                        for dot=1:obj.dots
                            subplot(floor(sqrt(obj.dots)), ceil(sqrt(obj.dots)), dot);
                            if (plot_st == 1)
                                obj.p(2*obj.dots+(dot-1)+1) = plot(t, Ps(dot));
                                ylim([0, 1]);
                            elseif (plot_st == 2)
                                obj.p(2*obj.dots+(dot-1)+1) = plot(t, 1-Pn(dot)+Ps(dot)+Pt(dot));
                                ylim([0, 2]);
                            end
                        end
                    end
                else
                    for dot=1:obj.dots
                        if (plot_style_lab == 1)
                            obj.p(2*(dot-1)+1).XData = [0, 2*Px(dot)-1];
                            obj.p(2*(dot-1)+1).YData = [0, 2*Py(dot)-1];
                            obj.p(2*(dot-1)+1).ZData = [0, 2*Pz(dot)-1];
                        elseif (plot_style_lab == 2)
                            obj.p(2*(dot-1)+1).XData = [obj.p(2*(dot-1)+1).XData, 2*Px(dot)-1];
                            obj.p(2*(dot-1)+1).YData = [obj.p(2*(dot-1)+1).YData, 2*Py(dot)-1];
                            obj.p(2*(dot-1)+1).ZData = [obj.p(2*(dot-1)+1).ZData, 2*Pz(dot)-1];
                        end
                        obj.p(2*(dot-1)+2).XData = [obj.p(2*(dot-1)+2).XData, 2*Pxr(dot)-1];
                        obj.p(2*(dot-1)+2).YData = [obj.p(2*(dot-1)+2).YData, 2*Pyr(dot)-1];
                        obj.p(2*(dot-1)+2).ZData = [obj.p(2*(dot-1)+2).ZData, 2*Pzr(dot)-1];
                        if (plot_st == 1)
                            obj.p(2*obj.dots+(dot-1)+1).XData = [obj.p(2*obj.dots+(dot-1)+1).XData, t];
                            obj.p(2*obj.dots+(dot-1)+1).YData = [obj.p(2*obj.dots+(dot-1)+1).YData, Ps(dot)];
                        elseif (plot_st == 2)
                            obj.p(2*obj.dots+(dot-1)+1).XData = [obj.p(2*obj.dots+(dot-1)+1).XData, t];
                            obj.p(2*obj.dots+(dot-1)+1).YData = [obj.p(2*obj.dots+(dot-1)+1).YData, 1-Pn(dot)+Ps(dot)+Pt(dot)];
                        end
                    end
                end
                
            end
            
            % Hamiltonian Property Getters/Setters
            % set(value)        : set to all dots (could be 1 dot)
            % set(dot, value)   : set to single dot
            % get()             : get value (for 1 dot, so dot 1)
            % get(dot)          : get value of single dot
            function setLarmorFrequency(obj, varargin)
                if (nargin == 2)
                    obj.larmor_frequency = varargin{1} * ones(1, obj.dots);
                else
                    dot = varargin{1};
                    obj.larmor_frequency(dot) = varargin{2};
                end
            end

            function larmor_frequency = getLarmorFrequency(obj, varargin)
                if (nargin == 1)
                    larmor_frequency = obj.larmor_frequency(1);
                else
                    dot = varargin{1};
                    larmor_frequency = obj.larmor_frequency(dot);
                end
            end

            function setRabiFrequency(obj, varargin)
                if (nargin == 2)
                    obj.rabi_frequency = varargin{1} * ones(1, obj.dots);
                else
                    dot = varargin{1};
                    obj.rabi_frequency(dot) = varargin{2};
                end
            end

            function rabi_frequency = getRabiFrequency(obj, varargin)
                if (nargin == 1)
                    rabi_frequency = obj.rabi_frequency(1);
                else
                    dot = varargin{1};
                    rabi_frequency = obj.rabi_frequency(dot);
                end
            end

            function setChargingEnergy(obj, varargin)
                if (nargin == 2)
                    obj.charging_energy = varargin{1} * ones(1, obj.dots);
                else
                    dot = varargin{1};
                    obj.charging_energy(dot) = varargin{2};
                end
            end

            function charging_energy = getChargingEnergy(obj, varargin)
                if (nargin == 1)
                    charging_energy = obj.charging_energy(1);
                else
                    dot = varargin{1};
                    charging_energy = obj.charging_energy(dot);
                end
            end
            
            function setSingletTripletEnergy(obj, varargin)
                if (nargin == 2)
                    obj.singlet_triplet_energy = varargin{1} * ones(1, obj.dots);
                else
                    dot = varargin{1};
                    obj.singlet_triplet_energy(dot) = varargin{2};
                end
            end

            function singlet_triplet_energy = getSingletTripletEnergy(obj, varargin)
                if (nargin == 1)
                    singlet_triplet_energy = obj.singlet_triplet_energy(1);
                else
                    dot = varargin{1};
                    singlet_triplet_energy = obj.singlet_triplet_energy(dot);
                end
            end
            
            % Hamiltonian Control Setters
            function setMicrowaveControl(obj, varargin)
                if (nargin == 2)
                    obj.microwave_control = varargin{1} * ones(1, obj.dots);
                else
                    dot = varargin{1};
                    obj.microwave_control(dot) = varargin{2};
                end
            end

            function setDetuningControl(obj, varargin)
                if (nargin == 2)
                    obj.detuning_control = varargin{1} * ones(1, obj.dots);
                else
                    dot = varargin{1};
                    obj.detuning_control(dot) = varargin{2};
                end
            end

            function setTunnelControl(obj, varargin)
                if (nargin == 2)
                    obj.tunnel_control = varargin{1} * ones(obj.dots, obj.dots);
                else
                    dota = varargin{1};
                    dotb = varargin{2};
                    obj.tunnel_control(dota, dotb) = varargin{3};
                    obj.tunnel_control(dotb, dota) = varargin{3};
                end
            end
            
        end
    
end
