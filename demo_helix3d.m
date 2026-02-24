function demo_helix3d()
    % -------- Demo parameters --------
    p0     = [0, 0, 0];      % start point
    p1     = [0, 0, 1.5];    % end point
    turns  = 6;              % number of turns
    r0     = 0.02;           % start radius
    r1     = 5;           % end radius (linearly changes to this)
    npts   = 3000;           % sampling points along the curve

    % -------- Generate ----------
    [X,Y,Z] = helix3d(p0, p1, turns, r0, r1, npts);

    % -------- Plot --------------
    figure; hold on; axis equal; grid on;
    plot3(X, Y, Z, 'LineWidth', 1.5);
    plot3([p0(1) p1(1)], [p0(2) p1(2)], [p0(3) p1(3)], 'k--', 'LineWidth', 1); % axis line
    xlabel('X'); ylabel('Y'); zlabel('Z'); title('3D Helix with Variable Radius');
    view(35, 25);
end

function [X,Y,Z] = helix3d(p0, p1, turns, r0, r1, npts)
    % helix3d
    % Inputs:
    %   p0, p1 : 1x3 vectors, start and end points (define helix axis)
    %   turns  : scalar, number of turns
    %   r0,r1  : start/end radius (linear interpolation along axis)
    %   npts   : number of samples along the curve
    % Outputs:
    %   X,Y,Z  : 1xn arrays (curve coordinates)

    % Axis direction and length
    d = p1 - p0;
    L = norm(d);
    if L < 1e-12
        error('Start and end points are too close.');
    end
    u = d / L;  % unit axis

    % Build an orthonormal basis (v1, v2) orthogonal to u
    % Choose a helper vector not parallel to u
    if abs(dot(u, [0 0 1])) < 0.99
        a = [0 0 1];
    else
        a = [0 1 0];
    end
    v1 = cross(u, a); v1 = v1 / norm(v1);
    v2 = cross(u, v1); % already unit-length if u and v1 are

    % Parameter t in [0,1]
    t = linspace(0, 1, npts);

    % Angle and radius along the helix
    theta  = 2*pi*turns*t;           % total turns
    radius = r0 + (r1 - r0)*t;       % linear radius change

    % Centerline along u
    C = p0 + (L*t').*u;              % n x 3

    % Circular offset in the orthogonal plane
    offs = (radius.*cos(theta))'.*v1 + (radius.*sin(theta))'.*v2;  % n x 3

    % Final coordinates
    P = C + offs;
    X = P(:,1).'; Y = P(:,2).'; Z = P(:,3).';
end