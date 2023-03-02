function [F] = maxwell(F_old, k, d, x_dot, time_step)

F = F_old + time_step*( k*x_dot - k*F_old/d );

% F = rk4(@(F_old)( k*x_dot - k*F_old/d ) , F_old, time_step);




end