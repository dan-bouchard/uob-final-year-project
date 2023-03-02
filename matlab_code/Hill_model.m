function [F_M, F_elements] = Hill_model(l_CE, l_M, l_CE_dot)

% Hill Model

CE.Fmax = 20;
CE.l_CEopt = 0.1;
CE.width = 0.56;
CE.f_vmax = 1.6;
CE.k_CE1 = 0.25;
CE.k_CE2 = 0.06;
CE.tau_c = 0.1; % time constant

PEE.l_PEE0 = 0.9*CE.l_CEopt;
PEE.F_PEE = 1.5;
%PEE.k_PEE
SDE.d_SDE = 0.3;
SEE.l_SEE0 = 0.2;
SEE.U0 = 0.04;

% Contractile element force
% Isometric force (Force length relation)
c = -1/(CE.width.^2);
CE_l_rel = l_CE/CE.l_CEopt;

if l_CE <= CE.l_CEopt*(1 - CE.width)
    F_isom = 0;
elseif l_CE >= CE.l_CEopt*(1 + CE.width)
    F_isom = 0;
else
    F_isom = c*CE_l_rel.^2 - 2*c*CE_l_rel + c + 1; 
end

 
% F_isom = CE.Fmax*f_L;

% Force velocity relation
v_max = CE.l_CEopt/CE.tau_c;
CE_l_dot_rel = l_CE_dot/v_max;

% concentric contraction
if CE_l_dot_rel <= -1
    f_v = 0;
elseif CE_l_dot_rel > -1 && CE_l_dot_rel <= 0
    f_v = (1 + CE_l_dot_rel)./(1 - CE_l_dot_rel/CE.k_CE1);
% eccentric contraction
else
    f_v = (1 + CE_l_dot_rel*(CE.f_vmax/CE.k_CE2))./(1 + CE_l_dot_rel/CE.k_CE2);
end    

F_CE = CE.Fmax*F_isom*f_v;

% Force of the parallel elastic element PEE
k_PEE = PEE.F_PEE*(CE.Fmax./(CE.l_CEopt*(CE.width + 1 - 0.9)).^3);

if l_CE >= PEE.l_PEE0
    F_PEE = k_PEE*(l_CE-PEE.l_PEE0).^3;
else % shorter than slack length
    F_PEE = 0;
end

% Force of the series elastic element SEE
k_SEE = CE.Fmax/(SEE.U0*SEE.l_SEE0).^3;
l_SEE = abs(l_M-l_CE);

if l_SEE >= SEE.l_SEE0
    F_SEE = k_SEE*(l_SEE-SEE.l_SEE0).^3;
else % shorter than slack length
    F_SEE = 0;
end


% Force of the serial damping element
F_SDE = SDE.d_SDE*l_CE_dot;

F_M = F_SEE;

F_elements = [F_CE F_isom f_v F_PEE F_SEE F_SDE]';


end