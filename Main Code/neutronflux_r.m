function r_prime = neutronflux_r(b,theta)
    syms x y
    m1 = tan(deg2rad(29.9999));
    eqn1 = y == (m1*x); %Constant boundary line    
    m2 = -tan(theta);
    eqn2 = y == (m2*x) + b;
    sol = solve([eqn1, eqn2], [x, y]);
    X1 = vpa(sol.x);
    Y1 = vpa(sol.y);
    X2 = 206.35;
    Y2 = (m2*X2) + b;
    DistanceArray = double([X1,Y1;X2,Y2]);
    r_prime = pdist(DistanceArray,'euclidean');
end