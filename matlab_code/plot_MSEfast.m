function plot_MSEfast(N)

% N = [2 5 10 15 25 40 60 100];
Errors = zeros(1, length(N));

for i=1:length(N)
    load(['files/learningTask2Fast/Error_N_', num2str(N(i))]);
    Errors(i)=Error;
end

figure;
plot(N,Errors,'x');
xlabel('N'); ylabel('MSE'); title('Mean-squared errors for different N values');
% set(gca, 'XScale', 'log');