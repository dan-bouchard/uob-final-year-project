function [net_after,sim_data] = simulate_hysteresis_sys (net,input,output)
% simulating mass-spring networks as used in [1] and [2]
%
% if output is given too, assumes teacher forcing is wanted
% Note: reads from net.init_dat.readout_type what readout to choose
%  
%  input: net		net_structure - has to be initialized by init_real_sd_net
%         input 	num_of_timesteps x num_of_input  input matrix to the system
%  		  output	used for teacher forcing 
%  output:
%  		 net_after	net structure after simulation (e.g. when m-file used at every single step)
%  		 sim_data  	data structure with all data harvested during simulation
%  
% 
%  helmut.hauser@bristol.ac.uk
%
%  [1] Hauser, H.; Ijspeert, A.; Füchslin, R.; Pfeifer, R. & Maass, W.
%  "Towards a theoretical foundation for morphological computation with compliant bodies"
%  Biological Cybernetics, Springer Berlin / Heidelberg, 2011, 105, 355-370 
%  http://www.springerlink.com/content/j236312507300638/
% 
%  [2] Hauser, H.; Ijspeert, A.; Füchslin, R.; Pfeifer, R. & Maass, W.
%  "The role of feedback in morphological computation with compliant bodies"
%  Biological Cybernetics, Springer Berlin / Heidelberg, 2012, 106, 595-613
%  http://www.springerlink.com/content/d54t39qh28561271/





% indices of input nodes
in_idx = find(net.W_in~=0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting necessary out net data structure
P = net.P;
W = net.W;
%  output_idx = net.output_idx;  %
time_step = net.init_data.time_step;
show_steps = net.init_data.show_steps;
sim_time = size(input,1)*time_step;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num = net.init_data.num;  % for the size of the data matrices
len = size(input,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (net.init_data.save_sim_data==1)
	% data matrices over time all together in sim_data structure
	sim_data.Fx = zeros(len,num);
	sim_data.Fy = zeros(len,num);

	sim_data.Sx_off = zeros(len,num);  % minus the offset
	sim_data.Sy = zeros(len,num);
	sim_data.Sxd = zeros(len,num);
	sim_data.Syd = zeros(len,num);
    sim_data.Sz = zeros(len,size(W.k1,1));
    sim_data.Sext = zeros(len,size(W.k1,1));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% now init values with values from the net
	sim_data.Sx(1,:) = P.states(:,1)';	% positions
	sim_data.Sy(1,:) = P.states(:,2)';
	sim_data.Sxd(1,:) = P.states(:,3)';	% velocities
	sim_data.Syd(1,:) = P.states(:,4)';
end
	% this data is put out in any way
	sim_data.O = zeros(len,net.init_data.nOutputs); % output weights * lengths of every spring
	% internal state - either D or Sx depending on "readout_type"
	sim_data.D  = zeros(len,size(W.k1,1)); % lengths of every spring
    
	sim_data.Sx = zeros(len,num); % x-positions of every mass point
    
%  	sim_data.O_off = zeros(len,net.init_data.nOutputs); % without offset 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulation loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx=0;  
for i=1:len
 	idx=idx+1; 	
 	
 	if(mod(idx,show_steps)==0) % when to show steps
 	  disp([' i = ' ,num2str(idx) , ' of ' , num2str(sim_time/time_step) ]);
 	end
 	
 	% set all old forces to zero (to get no unwanted acculumation)
 	P.force(:,1:2) = zeros(num,2);
 	
 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   	% go through all connections and calculate force
   	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for c=1:size(W.k1,1) % number of springs
	   % get actual points which are connected by this spring    
       from = W.from(c,1);
       to   = W.to(c,1);
       % actual distance with normed direction
       p_from = [ P.states(from,1) , P.states(from,2) ]';
       p_to   = [ P.states(to,1) , P.states(to,2) ]';
       [d,ndir] = e_distance( p_from,p_to); % distance and norm direction
	   %adding noise to the distance
%  	   d = d + net.dist_noise*rand(1,1);%-net.dist_noise*0.5;
	   if (net.init_data.save_sim_data==1 | strcmp(net.readout_type,'LENGTHS'))
			sim_data.D(idx,c)= d + net.dist_noise*rand(1,1); % sim_data.D is lengths
       end

	   
	% force amplitudes hysteresis
       x = d-W.l0(c,1); % extension
       x_dot = (d-W.dist_old(c,1))/time_step; % d/dt extension
	   [F, z_new] = Bouc_Wen_hysteresis(time_step, x, x_dot, W.z_old(c,1), W.beta(c,1), W.gamma(c,1), W.k3(c,1), W.k1(c,1), W.d3(c,1), W.d1(c,1));
       % [F, x_2, y_2, z_2] = Bouc_Wen_hysteresis(time_step, x, y, W.z_old(c,1) , W.beta(c,1), W.gamma(c,1));
	    
	  
        % add forces to the mass points which are not fixed
	   if(P.fixed(to,1)==0)   
	   	P.force(to,1)   = P.force(to,1) + (-1)*F*ndir(1,1); % f_x
	   	P.force(to,2)   = P.force(to,2) + (-1)*F*ndir(2,1); % f_y
	   end
	   if(P.fixed(from,1)==0)   
		P.force(from,1) = P.force(from,1) + (+1)*F*ndir(1,1); % f_x
	    P.force(from,2) = P.force(from,2) + (+1)*F*ndir(2,1); % f_y
	   end

	   % forces data
	   if (net.init_data.save_sim_data==1)
	   	sim_data.Fx(idx,to) = P.force(to,1);
	   	sim_data.Fy(idx,to) = P.force(to,2);
	   	sim_data.Fx(idx,from) = P.force(from,1);
	   	sim_data.Fy(idx,from) = P.force(from,2);
        sim_data.Sz(idx,c)    = W.z_old(c,1);
        sim_data.Sext(idx,c)  = x;
	   end
	   % update old distance with actual distnace
	   W.dist_old(c,1) = d; % update old distance of each spring
       W.z_old(c,1) = z_new;
    end 
    	   
    	   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add input signals
    % right now as Fx !! and on dimensional (exception => symmetric case!)
    % horizontal forcing
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    if (idx==1) % to make sure to get a nonzero idx	
		P.force(:,1) = P.force(:,1) + net.W_in * input(idx,:)' + net.W_fb*sim_data.O(idx,:)'; % x-dimension
			
    else
    	P.force(:,1) = P.force(:,1) + net.W_in * input(idx,:)' + net.W_fb*sim_data.O(idx-1,:)'; % x-dimension
    
    end
 
    	   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	   
    % get rid of all velocities and forces for fixed points	   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	   
    P.states(net.fixed_idx,3:4) = zeros(length(net.fixed_idx),2);
    P.force(net.fixed_idx,1:2)  = zeros(length(net.fixed_idx),2);
  
    
    
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simulation all dynamical system (mass points)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	for i=1:size(P.states,1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	   	%%%%%%%%%%%%%% EULER integrations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	   	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		
	   	x = ode_simple_ms_sys(time_step,P.states(i,1:4),P.force(i,1:2));
	   	P.states(i,:) = x;
	   		   	
	   	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	   	 if (net.init_data.save_sim_data==1)
	   		% getting data 
	   		sim_data.Sx(idx,i)  = P.states(i,1) + net.pos_noise*rand(size(P.states(i,1)));
	   		sim_data.Sy(idx,i)  = P.states(i,2) + net.pos_noise*rand(size(P.states(i,2)));
	   		sim_data.Sxd(idx,i) = P.states(i,3);
	   		sim_data.Syd(idx,i) = P.states(i,4);
	  	elseif (strcmp(net.readout_type,'POSITIONS'))
	  		sim_data.Sx(idx,i) = P.states(i,1)+ net.pos_noise*rand(size(P.states(i,1)));
	   	end


	end
	
		% getting outputs (depending on the readout scheme)
        switch net.readout_type
			
			case 'POSITIONS'
        		sim_data.O(idx,:) = net.W_out' * P.states(:,1); %  assuming just to use x-dimensions        	
        	case 'LENGTHS'
        		sim_data.O(idx,:) = net.W_out' * sim_data.D(idx,:)'; % single scalar output, (weights*lengths of every spring)	
        		
        end % switch
    if(length(find(isnan(sim_data.O(idx,:)))>0))
		sim_data.ERROR = 'NaN - ERROR / unstable simulation';
		disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
		disp('NaN - ERROR / unstable simulation');
		disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
		net_after = 0;
   		return
    end
end 	


% updating states for the net to be sent back
net_after = net;
net_after.P = P; % dynamic information about the points
net_after.W = W; % dynamic information about the connections
