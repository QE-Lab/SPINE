%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dU = solver_analytical_xz(dimension, dtH)

	% Note: analytically solved, assumes the following form for dtH
	% dtH = [a,  b
	%        b, -a]
    assert(dimension == 2, 'The Hamiltonian is not suitable for this solver!');
    assert(isreal(dtH), 'The Hamiltonian is not suitable for this solver!');
	a = dtH(1, 1);
	b = dtH(1, 2);
    
    % dU = expm(-1i * dtH)
    tmp = sqrt(a*a + b*b);
    tmp_cos = cos(tmp);
    tmp_sin = -1i / tmp * sin(tmp);
    dU = [tmp_cos + a * tmp_sin, b * tmp_sin;
          b * tmp_sin, tmp_cos - a * tmp_sin];  
    
end
