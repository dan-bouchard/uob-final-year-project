function [F_v] = f_v(l_dot_rel,k_CE1,k_CE2,v_max)
% F_v \in [0,1]
% The force-velocity relationship

% concentric contraction
if l_dot_rel <= -1
    F_v = 0;
elseif l_dot_rel > -1 && l_dot_rel <= 0
    F_v = (1 + l_dot_rel)./(1 - l_dot_rel/k_CE1);
% eccentric contraction
else
    F_v = (1 + l_dot_rel*(v_max/k_CE2))./(1 + l_dot_rel/k_CE2);
end
end




