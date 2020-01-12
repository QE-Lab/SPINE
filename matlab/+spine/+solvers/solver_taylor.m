%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dU = solver_taylor(dimension, dtH)

    % Default accuracy settings that can be set
    global solver_taylor_accuracy;
    if (isempty(solver_taylor_accuracy))
        solver_taylor_accuracy = 5;
    end
    global solver_taylor_scaling;
    if (isempty(solver_taylor_scaling))
        solver_taylor_scaling = 7;
    end
    
    % Scale the matrix for more accurate calculations
    M_small = (-1i * dtH) / 2^solver_taylor_scaling;
    M_power = M_small;
    
    % Exponentiation using series expansion
    result = complex(eye(dimension));
    factorial_i = 1;
    for i=1:(solver_taylor_accuracy - 1)
        factorial_i = factorial_i * i;
        result = result + M_power / factorial(i);
        M_power = M_power * M_small;
    end
    
    % Scale back, using squaring
    for i=1:solver_taylor_scaling
        result = result * result;
    end
    dU = result;
    
end
