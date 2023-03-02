function [params, out] = CompareModels(input_params)
% compare models with step input

if nargin == 1 % number of function input arguments
    disp('Using input parameters');
    ms = input_params.ms;
    hys = input_params.hys;
    maxwell = input_params.maxwell;
    kelvin = input_params.kelvin;
    hill = input_params.hill;
else
    disp('New random parameters');
    
    %%%%%%% MS System %%%%%%%
    % ms.k1 = rand_in_range_exp([1 200],1);
    % ms.k3 = rand_in_range_exp([1 100],1);
    % ms.d1 = rand_in_range_exp([1 200],1);
    % ms.d3 = rand_in_range_exp([1 100],1);

    % ms.k1 = rand_in_range_exp([10 50],1);
    % ms.k3 = rand_in_range_exp([5 25],1);
    % ms.d1 = rand_in_range_exp([10 50],1);
    % ms.d3 = rand_in_range_exp([5 25],1);
    % 
    % ms.k1 = rand_in_range_exp([10 100],1);
    % ms.k3 = rand_in_range_exp([4 60],1);
    % ms.d1 = rand_in_range_exp([10 100],1);
    % ms.d3 = rand_in_range_exp([4 60],1);

    ms.k1 = rand_in_range_exp([10 100],1);
    ms.k3 = rand_in_range_exp([1 100],1);
%     ms.d1 = rand_in_range_exp([10 100],1);
%     ms.d3 = rand_in_range_exp([1 100],1);
    ms.d1 = rand_in_range_exp([1 10],1);
    ms.d3 = rand_in_range_exp([1 10],1);
    
%     ms.k1 = 20;
%     ms.k3 = 10;
%     ms.d1 = 15;
%     ms.d3 = 5;

    % ms.k1 = rand_in_range_exp([0.1 1],1);
    % ms.k3 = rand_in_range_exp([0.1 1],1);
    % ms.d1 = rand_in_range_exp([0.1 1],1);
    % ms.d3 = rand_in_range_exp([0.1 1],1);
    % 
    % ms.k1 = 1;
    % ms.k3 = 0.6;
    % ms.d1 = 1;
    % ms.d3 = 0.6;


    %%%%%%% HYS System %%%%%%%
    % hys.k1 = rand_in_range_exp([1 200],1);
    % hys.k3 = rand_in_range_exp([1 100],1);
    % hys.d1 = rand_in_range_exp([1 200],1);
    % hys.d3 = rand_in_range_exp([1 100],1);
    hys.k1 = ms.k1;
    hys.k3 = ms.k3;
    hys.d1 = ms.d1;
    hys.d3 = ms.d3;


    hys.beta = rand_in_range_exp([1 100],1);
    hys.beta = rand_in_range_exp([0.1 100],1).*(2*randi(2)-3);
    % hys.beta = 50;
    hys.gamma = rand_in_range([-5 5],1);
    hys.gamma = rand_in_range_exp([0.1 100],1).*(2*randi(2)-3);
    % hys.gamma = 50;
    hys.A = 1;
    hys.n = 1.2;
    hys.alpha = 0;


    %%%%%%% MAXWELL System %%%%%%%
    % maxwell.k1 = rand_in_range_exp([1 200],1);
    % maxwell.k3 = rand_in_range_exp([1 100],1);
    % maxwell.d1 = rand_in_range_exp([1 200],1);
    % maxwell.d3 = rand_in_range_exp([1 100],1);

    maxwell.k1 = ms.k1;
    maxwell.k3 = ms.k3;
    maxwell.d1 = ms.d1;
    maxwell.d3 = ms.d3;


    %%%%%%% KELVIN AorB System %%%%%%%
    % kelvin.ks1 = rand_in_range_exp([1 200],1);
    % kelvin.ks3 = rand_in_range_exp([1 100],1);
    % kelvin.kp1 = rand_in_range_exp([1 200],1);
    % kelvin.kp3 = rand_in_range_exp([1 100],1);
    % kelvin.d1 = rand_in_range_exp([1 200],1);
    % kelvin.d3 = rand_in_range_exp([1 100],1);

    
%     kelvin.ks3 = 0.5*ms.k3;
%     kelvin.kp1 = rand_in_range_exp([1 10],1);
%     kelvin.kp3 = rand_in_range_exp([1 10],1);
    kelvin.kp1 = rand_in_range_exp([10 100],1);
%     kelvin.kp1 = 7;
%     kelvin.kp3 = 0.5*rand_in_range_exp([1 100],1);
    kelvin.kp3 = rand_in_range_exp([1 100],1);
%     kelvin.kp3 = 3;
    
    kelvin.ks1 = ms.k1;
    kelvin.ks3 = ms.k3;
    kelvin.d1 = ms.d1;
    kelvin.d3 = ms.d3;
%     kelvin.d1 = rand_in_range_exp([1 10],1);
%     kelvin.d3 = rand_in_range_exp([1 10],1);


    %%%%%%% HILL System %%%%%%%
%     hill.CE.Fmax = 20;
    hill.CE.Fmax = 2;
%     hill.CE.l_CEopt = 0.3;
    hill.CE.width = 0.56;
    hill.CE.f_vmax = 1.5;
    hill.CE.k_CE1 = 0.25;
    hill.CE.k_CE2 = 0.06;
    hill.CE.tau_c = 0.1; % time constant
    hill.CE.L1 = 3; % resting length of CE, SDE, PEE
    
%     hill.ks1 = rand_in_range_exp([10 100],1);
%     hill.ks3 = rand_in_range_exp([1 100],1);
%     hill.d1 = rand_in_range_exp([1 10],1);
%     hill.d3 = rand_in_range_exp([1 10],1);
%     hill.kp1 = rand_in_range_exp([10 100],1);
%     hill.kp3 = rand_in_range_exp([1 100],1);

    hill.ks1 = kelvin.ks1;
    hill.ks3 = kelvin.ks3;
    hill.d1 = kelvin.d1;
    hill.d3 = kelvin.d3;
    hill.kp1 = kelvin.kp1;
    hill.kp3 = kelvin.kp3;
    
end
params.ms = ms;
params.hys = hys;
params.maxwell = maxwell;
params.kelvin = kelvin;
params.hill = hill;

save('saved_params','params');

time_step = 0.01;
end_time = 100;
t = linspace(0,end_time,end_time/time_step + 1);
p=2; B=1;
f_t = ones(size(t));
num_zeros = length(f_t)-(end_time/(2*time_step))+1;
num_ones = (end_time/(2*time_step));
f_t(num_ones:end) = zeros(1,num_zeros);

% f_t(2501:5000) = zeros(1,2500);
% f_t(5001:7500) = ones(1,2500);
% f_t(7501:10001) = zeros(1,2501);

% f_t = B*sin(p*t);
% df_t = B*p*cos(p*t);
f1=2.11; f2=3.73; f3=4.33;
% f_t = sin(2*pi*f1.*t).*sin(2*pi*f2.*t).*sin(2*pi*f3.*t);
u_t = sin(2*pi*f1.*t).*sin(2*pi*f2.*t).*sin(2*pi*f3.*t);
a = (u_t-min(u_t))./(max(u_t)-min(u_t));
A = -4;
%a_t = (exp(A*a)-1)/(exp(A)-1);
a_t = a;

MS.x = zeros(size(t));
MS.x_dot = zeros(size(t));
MS.F = zeros(size(t));
MS.Fd = zeros(size(t));
MS.Fs = zeros(size(t));

HYS.x = zeros(size(t));
HYS.x_dot = zeros(size(t));
HYS.z = zeros(size(t));
HYS.F = zeros(size(t));
HYS.x2 = zeros(size(t));
HYS.x1_dot = zeros(size(t));
HYS.x1 = zeros(size(t));

MAX.x = zeros(size(t));
MAX.x_dot = zeros(size(t));
MAX.x2 = zeros(size(t));
MAX.x2_dot = zeros(size(t));
MAX.F = zeros(size(t));
MAX.Fs = zeros(size(t));
MAX.Fd = zeros(size(t));

KELVINA.x = zeros(size(t));
KELVINA.x_dot = zeros(size(t));
KELVINA.x2 = zeros(size(t));
KELVINA.x2_dot = zeros(size(t));
KELVINA.F = zeros(size(t));
KELVINA.Fsee = zeros(size(t));
KELVINA.Fsde = zeros(size(t));
KELVINA.Fpee = zeros(size(t));

KELVINB.x = zeros(size(t));
KELVINB.x_dot = zeros(size(t));
KELVINB.x2 = zeros(size(t));
KELVINB.x2_dot = zeros(size(t));
KELVINB.F = zeros(size(t));
KELVINB.Fsee = zeros(size(t));
KELVINB.Fsde = zeros(size(t));
KELVINB.Fpee = zeros(size(t));

HILL.x = zeros(size(t));
HILL.x_dot = zeros(size(t));
HILL.x2 = zeros(size(t));
HILL.x2_dot = zeros(size(t));
HILL.F = zeros(size(t));
HILL.Fsee = zeros(size(t));
HILL.Fsde = zeros(size(t));
HILL.Fpee = zeros(size(t));
HILL.F_CE = zeros(size(t));
HILL.F_equil = zeros(size(t));

% MAX.F(1) = maxwell.k1*f_t(1)+maxwell.k3*f_t(1)^3;
% MAX.x_dot(1) = df_t(1);
% MAX.x(1) = f_t(1);
% % MAX.x2(1) = f_t(1);
for i=1:length(t)-1
    if(mod(i,2000)==0) % when to show steps
 	  disp([' i = ' ,num2str(i) , ' of ' , num2str(length(t)) ]);
    end
    
    MS.x_dot(i+1) = MS.x_dot(i) + time_step*(f_t(i)-MS.F(i));
    MS.x(i+1) = MS.x(i) + time_step*MS.x_dot(i);
    MS.F(i+1) = ms.k1*MS.x(i+1)+ms.k3*(MS.x(i+1))^3+ms.d1*MS.x_dot(i+1)+ms.d3*(MS.x_dot(i+1))^3;
%     MS.Fs(i+1) = ms.k1*MS.x(i+1)+ms.k3*(MS.x(i+1))^3;
%     MS.Fd(i+1) = ms.d1*MS.x_dot(i+1)+ms.d3*(MS.x_dot(i+1))^3;
    
    HYS.x_dot(i+1) = HYS.x_dot(i) + time_step*(f_t(i)-HYS.F(i));
    HYS.x(i+1) = HYS.x(i) + time_step*HYS.x_dot(i);
    HYS.z(i+1) = HYS.z(i) + time_step*(hys.A*HYS.x_dot(i)-hys.beta*abs(HYS.x_dot(i))*abs(HYS.z(i)).^(hys.n-1)*HYS.z(i)-hys.gamma*HYS.x_dot(i)*abs(HYS.z(i)).^hys.n);
%     HYS.F(i+1) = hys.d1*HYS.x_dot(i)+hys.d3*(HYS.x_dot(i))^3+hys.alpha*(hys.k1*HYS.x(i)+hys.k3*(HYS.x(i))^3)+(1-hys.alpha)*(k*z(i));
    HYS.F(i+1) = hys.d1*HYS.x_dot(i+1)+hys.d3*(HYS.x_dot(i+1))^3+hys.alpha*(hys.k1*HYS.x(i+1)+hys.k3*(HYS.x(i+1))^3)+(1-hys.alpha)*(hys.k1*HYS.z(i+1)+hys.k3*(HYS.z(i+1))^3);
    
    MAX.x_dot(i+1) = MAX.x_dot(i) + time_step*(f_t(i)-MAX.F(i));
    MAX.x(i+1) = MAX.x(i) + time_step*MAX.x_dot(i);
    
%     nonlinear model
    max_fun = @(X2_dot)(maxwell.d3*(MAX.x_dot(i) - X2_dot)^3 + maxwell.d1*(MAX.x_dot(i) - X2_dot) - maxwell.k3*MAX.x2(i)^3 - maxwell.k1*MAX.x2(i));
    MAX.x2_dot(i+1) = fzero(max_fun, MAX.x2_dot(i));
    MAX.x2(i+1) = MAX.x2(i) + time_step*MAX.x2_dot(i+1);
%     MAX.Fs(i+1) = maxwell.k3*MAX.x2(i+1)^3+maxwell.k1*MAX.x2(i+1);
%     MAX.Fd(i+1) = maxwell.d3*(MAX.x_dot(i+1) - MAX.x2_dot(i+1)).^3 + maxwell.d1*(MAX.x_dot(i+1) - MAX.x2_dot(i+1));
    MAX.F(i+1) = maxwell.k3*MAX.x2(i+1)^3+maxwell.k1*MAX.x2(i+1);

% %     % linear model for maxwell (just ignore)
% %     HYS.x_dot(i+1) = HYS.x_dot(i) + time_step*(f_t(i)-HYS.F(i));
% %     HYS.x(i+1) = HYS.x(i) + time_step*HYS.x_dot(i);
% % %     HYS.x2_dot(i+1) = (-maxwell.k1*HYS.x2(i))/maxwell.d1 + HYS.x_dot(i);
% % % %     HYS.x2(i+1) = HYS.x2(i) + time_step*HYS.x2_dot(i);
% % %     HYS.x1_dot(i+1) = maxwell.k1/maxwell.d1*(HYS.x(i)-HYS.x1(i));
% % %     HYS.x1(i+1) = HYS.x1(i) + time_step*HYS.x1_dot(i);
% % %     HYS.F(i+1) = maxwell.d1*(HYS.x1_dot(i+1));
% % %     %HYS.F(i) = maxwell.d1*(HYS.x_dot(i)-HYS.x2_dot(i));
% %     
    
%     MAX.F(i+1) = MAX.F(i) + time_step*(maxwell.k1*MAX.x_dot(i)-(maxwell.k1*MAX.F(i))/maxwell.d1);
%     MAX.x_dot(i+1) = MAX.x_dot(i) + time_step*(f_t(i)-MAX.F(i));
%     MAX.x(i+1) = MAX.x(i) + time_step*MAX.x_dot(i);
    
    KELVINA.x_dot(i+1) = KELVINA.x_dot(i) + time_step*(f_t(i)-KELVINA.F(i));
    KELVINA.x(i+1) = KELVINA.x(i) + time_step*KELVINA.x_dot(i);
    KELVINA.x2_dot(i+1) = KELVINA.x_dot(i) - (kelvin.ks1/kelvin.d1)*KELVINA.x2(i);
    %kelvina_fun = @(X2_dot)(kelvin.d3*(KELVINA.x_dot(i) - X2_dot)^3 + kelvin.d1*(KELVINA.x_dot(i) - X2_dot) - kelvin.ks3*KELVINA.x2(i)^3 - kelvin.ks1*KELVINA.x2(i));
%     kelvina_fun(KELVINA.x2_dot(i))
    %KELVINA.x2_dot(i+1) = fzero(kelvina_fun, KELVINA.x2_dot(i));
    KELVINA.x2(i+1) = KELVINA.x2(i) + time_step*KELVINA.x2_dot(i+1);
    % x1_dot_kelvina  = KELVINA.x_dot(i+1)-KELVINA.x2_dot(i+1);
    %Fs_kelvina = kelvin.ks1*KELVINA.x2(i+1)+kelvin.ks3*KELVINA.x2(i+1)^3;
    
    %KELVINA.Fsee(i+1) = Fs_kelvina;
    % Fd_kelvina = kelvin.d1*(x1_dot_kelvina) + kelvin.d3*(x1_dot_kelvina)^3;
    % KELVINA.Fsde(i+1) = Fd_kelvina;
    % F_see_kelvina = 0.5*(Fs_kelvina+Fd_kelvina); 
    %F_pee_kelvina = kelvin.kp1*KELVINA.x(i+1) + kelvin.kp3*KELVINA.x(i+1)^3; 
    %KELVINA.Fpee(i+1) = F_pee_kelvina;
    KELVINA.F(i+1) = kelvin.ks1*KELVINA.x2(i+1)+kelvin.kp1*KELVINA.x(i+1);
    
%     KELVINB.F(i+1) = KELVINB.F(i) + time_step*(-1/kelvin.ks1 +KELVINB.x_dot(i) + (kelvin.kp1*KELVINB.x(i))/kelvin.d1 - (KELVINB.F(i)/kelvin.d1)*(1+kelvin.kp1/kelvin.ks1));
%     
    KELVINB.x_dot(i+1) = KELVINB.x_dot(i) + time_step*(f_t(i)-KELVINB.F(i));
    KELVINB.x(i+1) = KELVINB.x(i) + time_step*KELVINB.x_dot(i);
    kelvinb_fun = @(X2_dot)(kelvin.d3*(KELVINB.x_dot(i) - X2_dot)^3 + kelvin.d1*(KELVINB.x_dot(i) - X2_dot) + kelvin.kp3*(KELVINB.x(i) - KELVINB.x2(i))^3 + kelvin.kp1*(KELVINB.x(i) - KELVINB.x2(i)) - kelvin.ks3*KELVINB.x2(i)^3 - kelvin.ks1*KELVINB.x2(i));
    KELVINB.x2_dot(i+1) = fzero(kelvinb_fun, KELVINB.x2_dot(i));
    KELVINB.x2(i+1) = KELVINB.x2(i) + time_step*KELVINB.x2_dot(i+1);
    x1_dot_kelvinb  = KELVINB.x_dot(i+1)-KELVINB.x2_dot(i+1);
    Fsee_kelvinb = kelvin.ks1*KELVINB.x2(i+1)+kelvin.ks3*KELVINB.x2(i+1)^3; KELVINB.Fsee(i+1) = Fsee_kelvinb;
    Fsde_kelvinb = kelvin.d1*(x1_dot_kelvinb) + kelvin.d3*(x1_dot_kelvinb)^3; KELVINB.Fsde(i+1) = Fsde_kelvinb;
    x1_kelvinb = KELVINB.x(i+1) - KELVINB.x2(i+1);
    Fpee_kelvinb = kelvin.kp1*x1_kelvinb + kelvin.kp3*x1_kelvinb.^3; KELVINB.Fpee(i+1) = Fpee_kelvinb;
%     KELVINB.F(i+1) = 0.5*(Fsee_kelvinb+Fsde_kelvinb+Fpee_kelvinb);
    KELVINB.F(i+1) = Fsee_kelvinb;

    % HILL MODEL %
% %     HILL.x_dot(i+1) = HILL.x_dot(i) + time_step*(-HILL.F(i));
% %     HILL.x(i+1) = HILL.x(i) + time_step*HILL.x_dot(i);
% %     A = hill.CE.Fmax*a_t(i)*
% %     f_L(l_CE,l_rel,c,opt, width);
% %     f_v(l_dot_rel, k_CE1, k_CE2, f_vmax);
% %     HILL.F(i+1) = HILL.F(i) + time_step*((hill.ks1/hill.d1)*(hill.kp1*HILL.x(i) + hill.d1*HILL.x_dot(i) - (1+hill.kp1/hill.ks1)*F(i) + A));
    
    
    HILL.x_dot(i+1) = HILL.x_dot(i) + time_step*(f_t(i)-HILL.F(i));
    HILL.x(i+1) = HILL.x(i) + time_step*HILL.x_dot(i);
    
    l_CE = hill.CE.L1 + (HILL.x(i)-HILL.x2(i));
    opt = hill.CE.L1;
    l_rel = l_CE/opt;
    c = -1/(hill.CE.width.^2);
    v_max = opt/hill.CE.tau_c;
    
    hill_func_handle = @(X2_dot) hill_force_equilib(l_CE, l_rel, c, opt, hill.CE.width, v_max, hill.CE.k_CE1, hill.CE.k_CE2, hill.CE.f_vmax, hill.CE.Fmax, a_t(i), hill.d1, hill.d3, hill.kp1, hill.kp3, hill.ks1, hill.ks3, HILL.x(i), HILL.x2(i), HILL.x_dot(i), X2_dot);
                                 
    HILL.x2_dot(i+1) = fzero(hill_func_handle, HILL.x2_dot(i));
    x2_dot = HILL.x2_dot(i+1);
    hill_force_equilib(l_CE, l_rel, c, opt, hill.CE.width, v_max, hill.CE.k_CE1, hill.CE.k_CE2, hill.CE.f_vmax, hill.CE.Fmax, a_t(i), hill.d1, hill.d3, hill.kp1, hill.kp3, hill.ks1, hill.ks3, HILL.x(i), HILL.x2(i), HILL.x_dot(i), HILL.x2_dot(i+1));
    
    HILL.x2(i+1) = HILL.x2(i) + time_step*HILL.x2_dot(i+1);
    %x1_dot_hill  = HILL.x_dot(i+1)-HILL.x2_dot(i+1);
%     HILL.Fsee(i+1) = hill.ks1*HILL.x2(i+1)+hill.ks3*HILL.x2(i+1)^3;
    HILL.Fsee(i+1) = SEE_force_hill(HILL.x2(i+1), hill.ks1, hill.ks3);
%     HILL.Fsde(i+1) = hill.d1*(x1_dot_hill) + hill.d3*(x1_dot_hill)^3;
    x1_hill = HILL.x(i+1) - HILL.x2(i+1);
%     HILL.Fpee(i+1) = hill.kp1*x1_hill + hill.kp3*x1_hill.^3;
    HILL.Fpee(i+1) = PEE_force_hill(x1_hill, hill.kp1, hill.kp3);
    
    l_CE = hill.CE.L1 + (HILL.x(i+1)-HILL.x2(i+1));
    l_rel = l_CE/opt;
    F_isom = f_L(l_CE, l_rel, c, opt, hill.CE.width);
    %l_dot_rel = (HILL.x_dot(i+1)-HILL.x2_dot(i+1))./v_max;
    %F_v = f_v(l_dot_rel, hill.CE.k_CE1, hill.CE.k_CE2, v_max);
    %HILL.F_CE(i+1) = hill.CE.Fmax*a_t(i+1)*F_isom*F_v;
    %HILL.F_equil(i+1) = hill_force_equilib(l_CE, l_rel, c, opt, hill.CE.width, v_max, hill.CE.k_CE1, hill.CE.k_CE2, hill.CE.f_vmax, hill.CE.Fmax, a_t(i+1), hill.d1, hill.d3, hill.kp1, hill.kp3, hill.ks1, hill.ks3, HILL.x(i+1), HILL.x2(i+1), HILL.x_dot(i+1), HILL.x2_dot(i+1));
    
    HILL.F(i+1) = HILL.Fsee(i+1);
    % HILL.F(i+1) = HILL.Fsde(i+1) + HILL.Fpee(i+1) + HILL.F_CE(i+1);
end

transient = 40;
% figure; plot(MS.x(transient/time_step:end),MS.F(transient/time_step:end));

% figure; subplot(2,1,1); plot(MS.x,MS.F); hold on;
% plot(HYS.x,HYS.F); plot(MAX.x, MAX.F); plot(KELVINA.x, KELVINA.F); plot(KELVINB.x, KELVINB.F); hold off;
% legend('ms','hys','maxwell','kelvin A','kelvin B'); xlabel('x'); ylabel('F'); xlim([0 1]);
% subplot(2,1,2); plot(MS.x,MS.F); hold on;
% plot(HYS.x,HYS.F); plot(MAX.x, MAX.F); plot(KELVINA.x, KELVINA.F); plot(KELVINB.x, KELVINB.F); hold off;
% xlabel('x'); ylabel('F'); xlim([0 max(KELVINB.x)]);

% figure; plot(MAX.x, MAX.F); hold on; plot(MS.x, MS.F); hold off; legend('maxwell','voigt');
% figure; plot(t,MS.F); hold on; plot(t,f_t); plot(t,HYS.F); hold off; 
% legend('Voigt','f(t)','force derivation - dx1_dot');
% figure; plot(MAX.Fs); hold on; plot(MAX.Fd); legend('Fs','Fd');


% figure; plot(HYS.x,HYS.z); xlabel('x'); ylabel('z');
% y = linspace(-1,1,200);
% f = ms.k1*y+ms.k3*y.^3;
% figure; plot(y,f);

out.hys = HYS;
out.ms = MS;
out.maxwell = MAX;
out.kelvin_a = KELVINA;
out.kelvin_b = KELVINB;
out.hill = HILL;

end