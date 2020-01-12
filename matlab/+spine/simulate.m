%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function simulate(dimension, inHamiltonian, outOperation, solver, varargin)
    if (nargin < 5)
        % Operation simulation, allocate space for the operation
        U = complex(eye(dimension));
    else
        % State simulation, state passed as argument
        U = varargin{1};
    end

	% Let the user set the Hamiltonian and timestep
	[run, H, timestep] = inHamiltonian();
    while run

        % Determine (dt * H), the argument of expm()
        dtH = timestep * H;
        
		% Solve for the incremental operation
		dU = solver(dimension, dtH);

		% Update the overall operation
        U = dU * U;

		% Callback with the resulting operation
		outOperation(U);
        
        % Let the user set the Hamiltonian and timestep
        [run, H, timestep] = inHamiltonian();
    end
    
end
