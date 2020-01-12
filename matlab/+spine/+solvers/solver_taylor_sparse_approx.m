%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dU = solver_taylor_sparse_approx(dimension, dtH)

    % Default accuracy settings that can be set
    global solver_taylor_accuracy;
    if (isempty(solver_taylor_accuracy))
        solver_taylor_accuracy = 5;
    end
    global solver_taylor_scaling;
    if (isempty(solver_taylor_scaling))
        solver_taylor_scaling = 7;
    end
    global solver_taylor_tolerance;
    if (isempty(solver_taylor_tolerance))
        solver_taylor_tolerance = 1e-9;
    end
    
    % Scale the matrix for more accurate calculations
    M_small = (-sparse(dtH)) / 2^solver_taylor_scaling;
    M_power = M_small;
    
    % Exponentiation using series expansion
    result = sparse(complex(eye(dimension)));
    coeff = 1i;
    for m=2:solver_taylor_accuracy
        result = result + M_power * coeff;
        M_power = M_power * M_small;
        coeff = coeff * 1i / m;
    end
    
    % Scale back, using squaring
    for m=1:solver_taylor_scaling
        result = sparse(result .* (abs(result).^2 > solver_taylor_tolerance / 2^(solver_taylor_scaling - m + 1)));
        result = result * result;
    end

    dU = result;
    
end
