function [F_SEE] = SEE_force_hill(ext, k1, k3)
if ext >= 0
    F_SEE = k1*ext + k3*(ext).^3;
else % shorter than slack length
    F_SEE = 0;
end
end