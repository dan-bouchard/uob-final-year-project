function [dat]=learningTask(N, num_examples)

% N=100; num_examples=500000


% 2pif*timestep, timestep =1ms
n = 1:num_examples;
t = n/1000;


f1=2.11; f2=3.73; f3=4.33;
u = sin(2*pi*f1.*t).*sin(2*pi*f2.*t).*sin(2*pi*f3.*t);


sum1 = u;
if N>0
    for i=1:N
        %u_N = circshift(u,i);
        %u_N(1:i)=u_N(i+1);
        %sum1 = sum1+u_N;
        u_N = circshift(u,i);
        u_N(1:i)=u_N(i+1);
        % u_N = u_N.^2;
        sum1 = sum1+u_N;
    end
end
y=sum1;
% Y = mapstd(Y');
y = y./(N+1);
y = mapstd(y);

U = (u+1)/2;
A = -4;
a = (exp(A*U)-1)/(exp(A)-1);
% figure; plot(a); hold on; plot(U); legend('a','u'); hold off;


dat.u = u';
dat.y = y';

% figure; plot(t,u); hold on; plot(t,y); legend('u','y');