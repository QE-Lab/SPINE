%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dU = solver_diagonalization(~, dtH)

    % NOTE: the matrix is assumed to be real and symmetric!
    assert(isreal(dtH), 'The Hamiltonian is not suitable for this solver!')

    % 1) diagonalize the matrix: matrix = V * D * inv(V)
    [V, D] = eig(dtH, 'vector');
    
    % 2) exponentiate, exp(-i*D)
    D = exp(-1i * D);
    
    % 3) multiply V * D * inv(V)
    result =  V * diag(D) * V';
    dU = result;
    
end
