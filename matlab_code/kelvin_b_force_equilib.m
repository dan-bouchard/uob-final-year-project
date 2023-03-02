function [F_total] = kelvin_b_force_equilib(d1, d3, kp1, kp3, ks1, ks3, x, x2, x_dot, x2_dot)
% force equilibria for the kelvin_b model to setup function handle to calculate
% X2_dot

F_sde = d3*(x_dot - x2_dot).^3 + d1*(x_dot - x2_dot);
ext_pee = x-x2;
F_pee = PEE_force_hill(ext_pee, kp1, kp3);
% F_pee = kp3*(x - x2).^3 + kp1*(x - x2);
ext_see = x2;
F_see = SEE_force_hill(ext_see, ks1, ks3);
% F_see = ks3*x2^3 - ks1*x2;

F_total = F_sde+F_pee-F_see;

end

