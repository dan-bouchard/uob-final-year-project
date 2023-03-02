% Bouc_Wen_Model

time_step = 0.01;
end_time = 200;
t = linspace(0,end_time,end_time/time_step + 1);

% f1=2.11; f2=3.73; f3=4.33;
% u = sin(2*pi*f1.*t).*sin(2*pi*f2.*t).*sin(2*pi*f3.*t);

% Vectors setup
x = zeros(size(t));
y = zeros(size(t)); % x dot
z = zeros(size(t));
F = zeros(size(t));
% v = zeros(size(t));

% parameters
p=1; A=1; B=1;
alpha = 0.5; n = 1.2;
% k = 5; d = 10;
% k=1; d=0.5;

% k3 = 40; k1 = 19; k3 = 3.6; d1 = 4; d3 = 1.05;
% k1 = 17.2; k3 = 13.8; d1=91.8; d3 = 17.5;
k1 = 1; k3 = 1; d1=0.5; d3 = 0.5;

% doesn't seem to work with beta<0
beta = 0.25; gamma = 0.75;
beta = 0.3; gamma = -0.9;

f_t = B*sin(p*t);
% f_t = ones(size(t));

% ICS
% x(1) = 0; y(1) = 0; z(1) = 0; % v(1) = p;

for i=1:length(t)-1
%     X=x(i); Y=y(i); Z=z(i); % V=v(i);
%     U=u(i);
%     
%     x(i+1) = rk4(@(X)( Y ), X, time_step);
%     % y(i+1) = rk4(@(Y)( -d*Y-alpha*k*X-(1-alpha)*k*Z+B*sin(V) ), Y, time_step);
%     y(i+1) = rk4(@(Y)( -d*Y-alpha*k*X-(1-alpha)*k*Z+U ), Y, time_step);
%     % y(i+1) = rk4(@(Y)( -d3*Y^3-d1*Y-alpha*(k3*X^3+k1*X)-(1-alpha)*(k3*Z^3+k1*Z)+U ), Y, time_step);
%     z(i+1) = rk4(@(Z)( A*Y-beta*abs(Y)*abs(Z).^(n-1)*Z-gamma*Y*abs(Z).^n ) , Z, time_step);
%     %v(i+1) = rk4(@(V)( p ) , V, time_step);

%     y(i+1) = y(i) + time_step*(-d*y(i)-alpha*k*x(i)-(1-alpha)*k*z(i)+f_t(i));
    y(i+1) = y(i) + time_step*(f_t(i)-F(i));
    x(i+1) = x(i) + time_step*y(i);
    z(i+1) = z(i) + time_step*(A*y(i)-beta*abs(y(i))*abs(z(i)).^(n-1)*z(i)-gamma*y(i)*abs(z(i)).^n);
%     F(i+1) = d*y(i+1)+alpha*k*x(i+1)+(1-alpha)*k*z(i+1);
    F(i+1) = d1*y(i+1)+d3*(y(i+1))^3+alpha*(k1*x(i+1)+k3*x(i+1)^3)+(1-alpha)*(k1*z(i+1)+k3*z(i+1)^3);
    
        
end

transient = 50;

%F = alpha*k*x+(1-alpha)*k.*z+d.*y;
figure;
plot(x(transient/time_step:end),z(transient/time_step:end));
%plot(x,F);
xlabel('x'); ylabel('z'); title(['Bouc-Wen Hysteresis Loop (\beta=' num2str(beta) ', \gamma=' num2str(gamma) ')']);
% figure; plot(x,z); ylabel('z'); xlabel('x'); 
% plot(t,F); hold on; plot(t,u); xlim([10 15]);
% figure; plot(x(transient/time_step:end),z(transient/time_step:end));

%plot(t(10000:15000),F(10000:15000));
% plot(y(transient/time_step:end),z(transient/time_step:end));



