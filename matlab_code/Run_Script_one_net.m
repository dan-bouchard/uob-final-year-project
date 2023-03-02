
function Run_Script_one_net(n)

N = [10 15 25 40 60 100];

calculate_results_net(5,n);
for i=1:length(N)
    run_Net(N(i),n);
    calculate_results_net(N(i),n);
end

end