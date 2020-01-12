%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = fidelity(dimension, U, varargin)

    % Formula below is only valid for 2x2 operations
    assert(dimension == 2, 'The fidelity formula can only be used for 2x2 operations!');
    
    % Pass either Uideal, or theta and phi (Uideal will be calculated here)
    if (nargin < 4)
        Uideal = varargin{1};
    else
        theta = varargin{1};
        phi = varargin{2};
        Uideal = [cos(theta / 2), -1i * (cos(phi) + 1i * sin(phi)) * sin(theta / 2);
                  -1i * (cos(phi) - 1i * sin(phi)) * sin(theta / 2), cos(theta / 2)];
    end

    % Return the process fidelity
    F = abs(trace(Uideal' * U))^2/4;
    
end
