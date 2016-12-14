function X = solving_lines3(theta)
    % Finds the intersection points of the neutron ray(at defined angle) with
    % the neutron ray at critical angle for waveguide 2.
    r = 20;
    L = 200;
    %Initial Lines
    syms x y
    m1 = -tan(deg2rad(58.66)); %Constant boundary line
    b = -m1*456.35;
    eqn2 = y == (m1*x) + b + 30; 
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