function [X,Y,Z]=growing_spheres_demo(r_min,r_max,n_layers,n_points,randomize)
    % --------------------------------------
    % Parameters (you can modify these)
    % --------------------------------------
    % r_min    = 0.1;   % Minimum radius
    % r_max    = 1.0;   % Maximum radius
    % n_layers = 10;    % Number of concentric layers
    % n_points = 200;   % Points per layer
    % randomize = false; % true: random distribution, false: grid-like

    % --------------------------------------
    % Generate point cloud
    % --------------------------------------
    figure; hold on; axis equal; grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Growing Concentric Spheres (Points)');
    view(35,25);
    cmap = turbo(n_layers); % color for each layer

    % Radii list
    radii = linspace(r_min, r_max, n_layers);

    for k = 1:n_layers
        r = radii(k);
        % Generate angular coordinates
        if randomize
            theta = acos(2*rand(1,n_points)-1); % polar angle (0~π)
            phi   = 2*pi*rand(1,n_points);      % azimuth (0~2π)
        else
            n_side = ceil(sqrt(n_points));
            theta = linspace(0, pi, n_side);
            phi   = linspace(0, 2*pi, n_side);
            [Theta, Phi] = meshgrid(theta, phi);
            theta = Theta(:);
            phi   = Phi(:);
        end

        % Convert to Cartesian
        X = r * sin(theta) .* cos(phi);
        Y = r * sin(theta) .* sin(phi);
        Z = r * cos(theta);

        % Plot this layer of points
        scatter3(X, Y, Z, 15, 'MarkerFaceColor', cmap(k,:), ...
                 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
    end

    % Center
    plot3(0,0,0,'k.','MarkerSize',20);
end