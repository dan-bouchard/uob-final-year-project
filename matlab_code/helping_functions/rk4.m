function x2 = rk4(f, x1, deltaT)

    k1 = deltaT * f(x1);
    k2 = deltaT * f(x1 + k1/2);
    k3 = deltaT * f(x1 + k2/2);
    k4 = deltaT * f(x1 + k3);

    x2 = x1 + (k1 + k2*2 + 2*k3 + k4)/6;
end