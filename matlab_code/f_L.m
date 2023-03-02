function [F_isom] = f_L(l_CE,l_rel,c,opt, width)
% F_isom \in [0,1]
% The isometric force from the length relationship

if l_CE <= opt*(1 - width)
    F_isom = 0;
elseif l_CE >= opt*(1 + width)
    F_isom = 0;
else
    F_isom = c*l_rel.^2 - 2*c*l_rel + c + 1;
end

end

