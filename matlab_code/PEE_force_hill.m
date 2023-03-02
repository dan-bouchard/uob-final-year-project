function [F_PEE] = PEE_force_hill(ext, k1, k3)

if ext>= 0
    F_PEE = k1*ext + k3*(ext).^3;
else
    F_PEE = 0;
end
end