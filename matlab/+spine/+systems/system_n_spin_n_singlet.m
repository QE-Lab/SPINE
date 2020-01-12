%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef system_n_spin_n_singlet < spine.systems.system_spin
    
    properties
        size;
        storage_lookup;
        storage_size;
        
        % Hamiltonian Locations
        charging_energy_i;
        larmor_frequency_ip;
        larmor_frequency_im;
        tunnel_couplings_ip;
        tunnel_couplings_im;
        microwave_control_i;
        detuning_control_ip;
        detuning_control_im;
        
        % State Vector Locations
        probability_i;
        
    end
    
    methods
        
        function obj = system_n_spin_n_singlet(n_dots)
            
            % Determine the valid codes
            obj = obj@spine.systems.system_spin(n_dots);
            obj.size = 4^obj.dots;
            obj.storage_lookup = zeros(1, obj.size);
            index = 1;
            for i=0:obj.size-1
                if (obj.isValid(i))
                    obj.storage_lookup(i + 1) = index;
                    index = index + 1;
                end
            end
            obj.storage_size = index - 1;
    
            % Determine where the values are to be places in the Hamiltonian
            obj.charging_energy_i = cell(1, obj.dots);
            obj.larmor_frequency_ip = cell(1, obj.dots);
            obj.larmor_frequency_im = cell(1, obj.dots);
            obj.tunnel_couplings_ip = cell(obj.dots, obj.dots);
            obj.tunnel_couplings_im = cell(obj.dots, obj.dots);
            obj.microwave_control_i = cell(1, obj.dots);
            obj.detuning_control_ip = cell(1, obj.dots);
            obj.detuning_control_im = cell(1, obj.dots);
            
            % Determine where the measurements should take place
            obj.probability_i = cell(obj.dots, 4);

            for dot=0:obj.dots-1
                is00 = uint32(0 * 4^dot);
                is01 = uint32(1 * 4^dot);
                is02 = uint32(2 * 4^dot);
                is03 = uint32(3 * 4^dot);
                mask = uint32(3 * 4^dot);
                nmask = bitcmp(mask, 'uint32');
                for index=0:obj.size-1
                    storage_index = obj.storage_lookup(index + 1);
                    if (storage_index ~= 0)

                        % Add +tunnel_rates.At(dota, dotb) where dota is 0x01 and dotb is 0x02 and turns to 0x00 and 0x03 (and opposite), others don't care
                        % Add -tunnel_rates.At(dota, dotb) where dota is 0x02 and dotb is 0x01 and turns to 0x00 and 0x03 (and opposite), others don't care
                        if (dot < obj.dots-1)
                            for dotb=dot+1:obj.dots-1
                                is00b = uint32(0 * 4^dotb);
                                is01b = uint32(1 * 4^dotb);
                                is02b = uint32(2 * 4^dotb);
                                is03b = uint32(3 * 4^dotb);
                                maskb = uint32(3 * 4^dotb);
                                nmaskb = bitcmp(maskb, 'uint32');
                                if (bitand(index, mask) == is01) && (bitand(index, maskb) == is02b)
                                    otherindex = bitor(bitor(bitand(bitand(index, nmask), nmaskb), is00), is03b);
                                    storage_otherindex = obj.storage_lookup(otherindex + 1);
                                    obj.tunnel_couplings_ip{dot+1, dotb+1} = [obj.tunnel_couplings_ip{dot+1, dotb+1}, [storage_index; storage_otherindex]];
                                    otherindex = bitor(bitor(bitand(bitand(index, nmask), nmaskb), is03), is00b);
                                    storage_otherindex = obj.storage_lookup(otherindex + 1);
                                    obj.tunnel_couplings_ip{dot+1, dotb+1} = [obj.tunnel_couplings_ip{dot+1, dotb+1}, [storage_index; storage_otherindex]];
                                else
                                    if (bitand(index, mask) == is02) && (bitand(index, maskb) == is01b)
                                        otherindex = bitor(bitor(bitand(bitand(index, nmask), nmaskb), is00), is03b);
                                        storage_otherindex = obj.storage_lookup(otherindex + 1);
                                        obj.tunnel_couplings_im{dot+1, dotb+1} = [obj.tunnel_couplings_im{dot+1, dotb+1}, [storage_index; storage_otherindex]];
                                        otherindex = bitor(bitor(bitand(bitand(index, nmask), nmaskb), is03), is00b);
                                        storage_otherindex = obj.storage_lookup(otherindex + 1);
                                        obj.tunnel_couplings_im{dot+1, dotb+1} = [obj.tunnel_couplings_im{dot+1, dotb+1}, [storage_index; storage_otherindex]];
                                    end
                                end
                            end
                        end

                        % Add +larmor_frequency[dot]/2 where this dot is 0x02, others don't care
                        % Add -larmor_frequency[dot]/2 where this dot is 0x01, others don't care
                        if (bitand(index, mask) == is01)
                            obj.larmor_frequency_im{dot+1} = [obj.larmor_frequency_im{dot+1}, storage_index];
                        else
                            if (bitand(index, mask) == is02)
                                obj.larmor_frequency_ip{dot+1} = [obj.larmor_frequency_ip{dot+1}, storage_index];
                            end
                        end

                        % Add charging_energy[dot] where this dot is 0x03, others don't care
                        if (bitand(index, mask) == is03)
                            obj.charging_energy_i{dot+1} = [obj.charging_energy_i{dot+1}, storage_index];
                        end

                        % Add microwave_control[dot] where this dot is 0x01 and turns to 0x02 (and vice versa), others don't care
                        if (bitand(index, mask) == is01)
                            otherindex = bitor(bitand(index, nmask), is02);
                            storage_otherindex = obj.storage_lookup(otherindex + 1);
                            obj.microwave_control_i{dot+1} = [obj.microwave_control_i{dot+1}, [storage_index; storage_otherindex]];
                        end

                        % Add -detuning_control[dot] where this dot is 0x03, others don't care
                        if (bitand(index, mask) == is03)
                            obj.detuning_control_im{dot+1} = [obj.detuning_control_im{dot+1}, storage_index];
                            for d=0:obj.dots-1
                                if (d ~= dot)
                                    obj.detuning_control_ip{d+1} = [obj.detuning_control_ip{d+1}, storage_index];
                                end
                            end
                        end

                        % Store locations required to determine the measurement probabilities
                        if (bitand(index, mask) == is00)
                            obj.probability_i{dot+1, 1} = [obj.probability_i{dot+1, 1}, storage_index];
                        else
                            if (bitand(index, mask) == is01)
                                otherindex = bitor(bitand(index, nmask), is02);
                                storage_otherindex = obj.storage_lookup(otherindex + 1);
                                obj.probability_i{dot+1, 2} = [obj.probability_i{dot+1, 2}, storage_index];
                                obj.probability_i{dot+1, 3} = [obj.probability_i{dot+1, 3}, storage_otherindex];
                            else
                                if (bitand(index, mask) == is03)
                                    obj.probability_i{dot+1, 4} = [obj.probability_i{dot+1, 4}, storage_index];
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function V = isValid(obj, index)
            if (index >= obj.size)
                V = false;
            else
                electrons = 0;
                for i=0:obj.dots-1
                    code = bitand(index, 3);
                    if ((code == 1) || (code == 2))
                        electrons = electrons + 1;
                    else
                        if (code == 3)
                            electrons = electrons + 2;
                        end
                    end
                    index = floor(index / 4);
                end
                V = (electrons == obj.dots);
            end
        end
        
        function index = getIndex(obj, states)
            index = 0;
            for dot=0:obj.dots-1
                index = index + uint32(states(dot+1) * 4^dot);
            end
            index = obj.storage_lookup(index + 1);
        end
        
        function index = getIndexMeasurement(obj, dot, state)
            index = [];
            if ((dot <= obj.dots) && (state < 4))
                index = obj.probability_i{dot, state+1};
            end
        end
        
        function dimension = getDimension(obj)
            dimension = obj.storage_size;
        end

        function H = updateHamiltonian(obj)

            % Clear the entire Hamiltonian
            H = zeros(obj.storage_size, obj.storage_size);

            % Add the charging energies
            for dot=1:obj.dots
                value = obj.charging_energy(dot);
                indices = obj.charging_energy_i{dot};
                for i=1:length(indices)
                    H(indices(i), indices(i)) = H(indices(i), indices(i)) + value;
                end
            end

            % Add the larmor frequencies
            for dot=1:obj.dots
                value = obj.larmor_frequency(dot);
                indices_p = obj.larmor_frequency_ip{dot};
                for i=1:length(indices_p)
                    H(indices_p(i), indices_p(i)) = H(indices_p(i), indices_p(i)) + value/2;
                end
                indices_m = obj.larmor_frequency_im{dot};
                for i=1:length(indices_m)
                    H(indices_m(i), indices_m(i)) = H(indices_m(i), indices_m(i)) - value/2;
                end
            end

            % Add the tunnel couplings
            for dot=1:obj.dots
                if (dot < obj.dots)
                    for dotb=dot+1:obj.dots
                        value = obj.tunnel_control(dot, dotb);
                        indices_p = obj.tunnel_couplings_ip{dot, dotb};
                        for i=1:length(indices_p)
                            H(indices_p(1, i), indices_p(2, i)) = H(indices_p(1, i), indices_p(2, i)) + value;
                            H(indices_p(2, i), indices_p(1, i)) = H(indices_p(2, i), indices_p(1, i)) + value;
                        end
                        indices_m = obj.tunnel_couplings_im{dot, dotb};
                        for i=1:length(indices_m)
                            H(indices_m(1, i), indices_m(2, i)) = H(indices_m(1, i), indices_m(2, i)) - value;
                            H(indices_m(2, i), indices_m(1, i)) = H(indices_m(2, i), indices_m(1, i)) - value;
                        end
                    end
                end
            end

            % Add the microwave control
            for dot=1:obj.dots
                value = obj.microwave_control(dot) * obj.rabi_frequency(dot);
                indices = obj.microwave_control_i{dot};
                for i=1:length(indices)
                    H(indices(1, i), indices(2, i)) = H(indices(1, i), indices(2, i)) + value;
                    H(indices(2, i), indices(1, i)) = H(indices(2, i), indices(1, i)) + value;
                end
            end

            % Add the detuning control
            for dot=1:obj.dots
                value = obj.detuning_control(dot);
                indices_p = obj.detuning_control_ip{dot};
                for i=1:length(indices_p)
                    H(indices_p(i), indices_p(i)) = H(indices_p(i), indices_p(i)) + value;
                end
                indices_m = obj.detuning_control_im{dot};
                for i=1:length(indices_m)
                    H(indices_m(i), indices_m(i)) = H(indices_m(i), indices_m(i)) - value;
                end
            end
        end
        
        function toString(obj, filename_H, filename_state)

            % Create the storage space for the text
            H = cell(obj.storage_size, obj.storage_size);
            for i=1:obj.storage_size
                for j=1:obj.storage_size
                    H{i, j} = '0';
                end
            end

            % Add the charging energies
            for dot=1:obj.dots
                indices = obj.charging_energy_i{dot};
                for i=1:length(indices)
                    H{indices(i), indices(i)} = sprintf('%s %s%d', H{indices(i), indices(i)}, '+Uc_', dot);
                end
            end

            % Add the larmor frequencies
            for dot=1:obj.dots
                indices_p = obj.larmor_frequency_ip{dot};
                for i=1:length(indices_p)
                    H{indices_p(i), indices_p(i)} = sprintf('%s %s%d', H{indices_p(i), indices_p(i)}, '+1/2w0_', dot);
                end
                indices_m = obj.larmor_frequency_im{dot};
                for i=1:length(indices_m)
                    H{indices_m(i), indices_m(i)} = sprintf('%s %s%d', H{indices_m(i), indices_m(i)}, '-1/2w0_', dot);
                end
            end

            % Add the tunnel couplings
            for dot=1:obj.dots
                if (dot < obj.dots)
                    for dotb=dot+1:obj.dots
                        indices_p = obj.tunnel_couplings_ip{dot, dotb};
                        for i=1:length(indices_p)
                            H{indices_p(1, i), indices_p(2, i)} = sprintf('%s %s%d%d', H{indices_p(1, i), indices_p(2, i)}, '+t0_', dot, dotb);
                            H{indices_p(2, i), indices_p(1, i)} = sprintf('%s %s%d%d', H{indices_p(2, i), indices_p(1, i)}, '+t0_', dot, dotb);
                        end
                        indices_m = obj.tunnel_couplings_im{dot, dotb};
                        for i=1:length(indices_m)
                            H{indices_m(1, i), indices_m(2, i)} = sprintf('%s %s%d%d', H{indices_m(1, i), indices_m(2, i)}, '-t0_', dot, dotb);
                            H{indices_m(2, i), indices_m(1, i)} = sprintf('%s %s%d%d', H{indices_m(2, i), indices_m(1, i)}, '-t0_', dot, dotb);
                        end
                    end
                end
            end

            % Add the microwave control
            for dot=1:obj.dots
                indices = obj.microwave_control_i{dot};
                for i=1:length(indices)
                    H{indices(1, i), indices(2, i)} = sprintf('%s %s%d', H{indices(1, i), indices(2, i)}, 'mw_', dot);
                    H{indices(2, i), indices(1, i)} = sprintf('%s %s%d', H{indices(2, i), indices(1, i)}, 'mw_', dot);
                end
            end

            % Add the detuning control
            for dot=1:obj.dots
                indices_p = obj.detuning_control_ip{dot};
                for i=1:length(indices_p)
                    H{indices_p(i), indices_p(i)} = sprintf('%s %s%d', H{indices_p(i), indices_p(i)}, '+e_', dot);
                end
                indices_m = obj.detuning_control_im{dot};
                for i=1:length(indices_m)
                    H{indices_m(i), indices_m(i)} = sprintf('%s %s%d', H{indices_m(i), indices_m(i)}, '-e_', dot);
                end
            end

            % Store to file
            fileID = fopen(filename_H, 'w');
            width = zeros(1, obj.storage_size);
            for i=1:obj.storage_size
                for j=1:obj.storage_size
                    len = length(H{i, j});
                    if (len > width(j))
                        width(j) = len;
                    end
                end
            end
            for i=1:obj.storage_size
                for j=1:obj.storage_size
                    fprintf(fileID, sprintf('%%%ds', width(j)), H{i, j});
                    fprintf(fileID, ', ');
                end
                fprintf(fileID, '\r\n');
            end
            fclose(fileID);

            % Determine the states
            fileID = fopen(filename_state, 'w');
            for index=0:obj.size-1
                index_run = index;
                str = '    ';
                electrons = 0;
                for i=0:obj.dots-1
                    code = bitand(index_run, 3);
                    if (code == 0)
                        str(i+1) = '-';
                    else
                        str(i+1) = '0' + code-1;
                    end
                    if ((code == 1) || (code == 2))
                        electrons = electrons + 1;
                    else
                        if (code == 3)
                            electrons = electrons + 2;
                        end
                    end
                    index_run = floor(index_run / 4);
                end
                if (electrons == obj.dots)
                    fprintf(fileID, '%s\r\n', str);
                end
            end
            fclose(fileID);

        end

    end
    
end
