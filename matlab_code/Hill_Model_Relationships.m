% Hill_Model_Relationships

CE.Fmax = 20;
CE.l_CEopt = 0.1;
CE.width = 0.56;
CE.f_vmax = 1.5;
CE.k_CE1 = 0.25;
CE.k_CE2 = 0.06;
CE.tau_c = 0.1; % time constant

PEE.L_PEE0 = 0.85;
PEE.l_PEE0 = PEE.L_PEE0*CE.l_CEopt;
PEE.F_PEE = 1.3;
%PEE.k_PEE
SDE.d_SDE = 0.3;
SEE.l_SEE0 = 0.2;
SEE.U0 = 0.04;


l_CE = linspace(0,0.3,200);
l_M = linspace(0.2,0.42,200);
l_SEE = (l_M-CE.l_CEopt);
SEE_l_rel = l_SEE/(CE.l_CEopt+SEE.l_SEE0);

c = -1/(CE.width.^2);
CE_l_rel = l_CE/CE.l_CEopt;

% k_PEE = PEE.F_PEE*(CE.Fmax./(CE.l_CEopt*(CE.width + 1 - PEE.L_PEE0)).^3);
k_PEE = PEE.F_PEE*(CE.Fmax./(CE.l_CEopt*(CE.width + 1.3*(1 - PEE.L_PEE0))).^3)

% k_PEE = 90000;

l_PEE_max = PEE.l_PEE0 + (PEE.F_PEE*CE.Fmax/k_PEE).^(1/3);
k_SEE = CE.Fmax/(SEE.U0*SEE.l_SEE0).^3

F_isom = zeros(size(l_CE));
F_PEE  = zeros(size(l_CE));
F_SEE  = zeros(size(l_SEE));
for i=1:length(F_isom)
    F_isom(i) = F_isometric(l_CE(i),CE_l_rel(i),c,CE.l_CEopt, CE.width);
    F_PEE(i)  = PEE_force(l_CE(i), PEE.l_PEE0, k_PEE, l_PEE_max, PEE.F_PEE*CE.Fmax);
    F_SEE(i) = SEE_force(l_SEE(i),SEE.l_SEE0,k_SEE);
    F_SEE(i) = SEE_force(l_CE(i),SEE.l_SEE0,k_SEE);
    
end
F_isom = CE.Fmax*F_isom;
F_tot = F_isom + F_PEE;

figure; plot(l_CE,F_isom, l_CE,F_PEE, l_CE,F_tot);
% plot(l_CE,F_tot);

% Force velocity relation
v_max = CE.l_CEopt/CE.tau_c;


l_CE_dot = linspace(-3,3,300);
f_v = zeros(size(l_CE_dot));
CE_l_dot_rel = zeros(size(l_CE_dot));
for i=1:length(f_v)
    CE_l_dot_rel(i) = l_CE_dot(i)/v_max;
    f_v(i) = Vel_Rel(CE_l_dot_rel(i),CE.k_CE1,CE.k_CE2,CE.f_vmax);
end
figure; plot(CE_l_dot_rel,f_v);


figure; plot(SEE_l_rel,F_SEE);
figure; plot(l_CE,F_SEE);
% l_SEE = abs(l_M-l_CE);


% SEE.d = F_M/3 if F_M<=F_Max
% else SEE.d = F_Max/3


% Go with the model of EMG data 
% a = (exp(A*u(t))-1)/(exp(A)-1)
% normalise input signal u by taking abs
% Can explain about the other model of neural inputs

A = -4;


function [F_isom] = F_isometric(l_CE,l_rel,c,opt, width)
    
    if l_CE <= opt*(1 - width)
        F_isom = 0;
        
    elseif l_CE >= opt*(1 + width)
        F_isom = 0;
    else
        F_isom = c*l_rel.^2 - 2*c*l_rel + c + 1;
    end

end

function [F_PEE] = PEE_force(l_CE, l_0, k,l_max,F_max)
if l_CE >= l_0 && l_CE <= l_max
    F_PEE = k*(l_CE-l_0).^3;
elseif l_CE > l_max
    F_PEE = F_max;
else % shorter than slack length
    F_PEE = 0;
end
end

function [f_v] = Vel_Rel(l_dot_rel,k_CE1,k_CE2,v_max)
% concentric contraction
if l_dot_rel <= -1
    f_v = 0;
elseif l_dot_rel > -1 && l_dot_rel <= 0
    f_v = (1 + l_dot_rel)./(1 - l_dot_rel/k_CE1);
% eccentric contraction
else
    f_v = (1 + l_dot_rel*(v_max/k_CE2))./(1 + l_dot_rel/k_CE2);
end
end

function [F_SEE] = SEE_force(l_SEE,slack,k)
if l_SEE >= slack
    F_SEE = k*(l_SEE-slack).^3;
else % shorter than slack length
    F_SEE = 0;
end
end