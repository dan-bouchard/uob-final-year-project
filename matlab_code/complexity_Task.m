function [dat]=complexity_Task(N, num_examples)
% N=2; num_examples=500000

n = 1:num_examples;
t = n/1000;


f1=2.11; f2=3.73; f3=4.33;
u = sin(2*pi*f1.*t).*sin(2*pi*f2.*t).*sin(2*pi*f3.*t);

coeffs = (2*rand(N+1,1)-1);
sum1 = coeffs(1)*ones(size(t));

if N>0
    for i=1:N
        sum1 = sum1+coeffs(i+1)*u.^i;
    end
end

y=sum1;
% Y = mapstd(Y');
y = y./(N+1);
y = mapstd(y);

dat.u = u';
dat.y = y';
dat.coeffs = coeffs;

figure; plot(u); hold on; plot(y); legend('u','y');