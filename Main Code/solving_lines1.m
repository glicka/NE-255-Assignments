function X = solving_lines1(theta)
    % Finds the intersection points of the neutron ray(at defined angle) with
    % the neutron ray at critical angle for waveguide 1.
    r = 56.35;
    L = 150;
    %Initial Lines
    syms x y
    m1 = tan(deg2rad(29.9999));
    eqn2 = y == (m1*x) + 30; %Constant boundary line
    m2 = tan(deg2rad(theta));
    eqn3 = y == m2*x;
    % Solving the lines to find the intersection points.
    sol = solve([eqn2, eqn3], [x, y]);
    X = vpa(sol.x);
    y = vpa(sol.y);
    % Intersection of the two lines cannot be outside our device.
    if X > (r+L)
        X = r+L;
    end
end