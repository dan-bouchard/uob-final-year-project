N = [5 10 15 25 40 60 100];
All_Errors = zeros(10,length(N));
x = ones(10,length(N));

for i=1:length(N)
    load(['files/smallNets_kelvin_b/Error_kelvin_b_N_', num2str(N(i))]);
    All_Errors(:,i)=Error;
    x(:,i) = x(:,i)*N(i);
end
All_errors = reshape(All_Errors,10*length(N),1);
% idx = find(All_errors>0);
x = reshape(x,10*length(N),1);
figure;
plot(x,All_errors,'x');
xlabel('N'); ylabel('MSE'); title('Mean-squared errors for memory task - Kelvin B');
% set(gca, 'YScale', 'log');
% set(gca, 'XScale', 'log');