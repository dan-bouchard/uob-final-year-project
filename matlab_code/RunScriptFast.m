%Run Script

N = [2 5 10 20 40 60 100];
%N = [40 60 100];

for i=1:length(N)
    disp([' j = ' ,num2str(i) , ' of  ', num2str(length(N))]);
    runOneNetwork(N(i),'spring');
    calculateResultsFast(N(i));
end

plot_MSEfast(N);
