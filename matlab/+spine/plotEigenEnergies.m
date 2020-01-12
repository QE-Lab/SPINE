%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotEigenEnergies(system, varargin)

    % Number of points
    points = 1001;
    if (nargin == 2)
        points = varargin{1};
    end

    % Only for a double dot spin system
    assert(isa(system, 'spine.systems.system_spin'), 'Only derivatives of spine.systems.system_spin are accepted!');
    assert(system.dots == 2, 'Only double dot systems are accepted!');
    
    % Sweep the relative detuning from -U to U
    % detuning_control = obj.detuning_control(2) - obj.detuning_control(1);
    detuning_old_2 = system.detuning_control(2);
    detuning_old_1 = system.detuning_control(1);
    charging_energy = max(system.charging_energy(1), system.charging_energy(2));
    system.detuning_control(1) = 0;
    detuning = linspace(-1.25*charging_energy, 1.25*charging_energy, points);
    energies = zeros(system.getDimension(), points);
    for i=1:length(detuning)
        system.detuning_control(2) = detuning(i);
        H = system.updateHamiltonian();
        energies(:, i) = eig(H);
    end
    system.detuning_control(2) = detuning_old_2;
    system.detuning_control(1) = detuning_old_1;

    % Plot
    for i=1:system.getDimension()
        plot(detuning, energies(i, :), 'k'); hold on;
    end
    xlabel('Detuning control');
    ylabel('Energy');
    
end
