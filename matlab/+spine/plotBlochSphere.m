%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J.P.G. van Dijk                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotBlochSphere()

    % Draw a simple Bloch sphere
    hold on;
    axis off;
    view(120, 25);
    camva(6);
    theta = linspace(0, pi, 50);
    for phi = 0:0.1*pi:2*pi
        plot3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta), ':', 'color', 0.3*[1 1 1]);
    end
    phi = linspace(0, 2*pi, 100);
    for theta = 0:0.1*pi:pi
        plot3(sin(theta)*cos(phi), sin(theta)*sin(phi), ones(1, 100)*cos(theta), ':', 'color', 0.3*[1 1 1]);
    end
    plot3([-1, 1], [0, 0], [0, 0], '-', 'color', 0*[1 1 1]);
    plot3([0, 0], [-1, 1], [0, 0], '-', 'color', 0*[1 1 1]);
    plot3([0, 0], [0, 0], [-1, 1], '-', 'color', 0*[1 1 1]);
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([-1 1]);
    text(0, 0, -1, '|1>', 'fontsize', 18);
    text(0, 0, 1, '|0>', 'fontsize', 18);
    
end
