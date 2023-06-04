clear;
clc;
count = 0;
% INPUT parameters
% Bus data
% Bus No		P,MW		Q,MVAR		V_mag,PU	V_phase		
Bus1=   [1		 0.0		0.0			1.0			 0.0	;
         0          0       0           0               0;
         0          0       0           0               0;         
	 	 4       2.0     	1.6 		0.991		-0.370  ;
	 	 5       3.0     	0.4 		0.9888		-0.544  ;
	 	 6       2.0     	-0.4 		0.986		-0.697  ;
	 	 7       1.5     	1.2 		0.985		-0.704  ;
	 	 8       4.0     	2.7 		0.979		-0.763  ;
	 	 9       5.0     	0.8 		0.971		-1.451  ;
	 	 10      1.0     	0.9 		0.977		-0.770  ;
	 	 11      0.6     	-0.5		0.971		-1.525  ;
	 	 12      4.5     	-1.7 		0.969		-1.836  ;
	 	 13      1.0     	0.9 		0.994		-0.332  ;
	 	 14      1.0     	-1.1 		0.995		-0.459  ;
	 	 15      1.0     	0.9 		0.992		-0.527  ;
	 	 16      2.1     	-0.8		0.991		-0.596  ];
		   
% Branch data
%		Line	Bus		Bus	PU Branch 	PU Branch	
%		no:		1	    2	resistance	reactance	

Line=    [ 1      0      1       0       0           ;
           0      0     0       0       0           ;
           0      0     0       0       0           ;
           0      0     0       0       0           ;
           0      0     0       0       0           ;
           0      0     0       0       0           ;
           0      0     0       0       0           ;
           0      0     0       0       0           ;
           0      0     0       0       0           ;
           0      0     0       0       0           ;
          11	 1		4 	 0.075   	0.1     	;
          12     4   	5 	 0.08    	0.11    	;
          13     4   	6 	 0.09    	0.18    	;
          14     6   	7 	 0.04    	0.04    	;
          15     5   	11 	 0.04    	0.04    	;
          16     1   	8 	 0.11    	0.11    	;
          17     8   	10	 0.11    	0.11    	;
          18     8   	9 	 0.08    	0.11    	;          
          19     9   	11	 0.11    	0.11    	;
          20     9   	12	 0.08    	0.11    	;
          21     10  	14	 0.04    	0.04    	;
          22     1   	13	 0.11    	0.11    	;
          23     13  	15	 0.08    	0.11    	;
          24     13  	14	 0.09    	0.12    	;          
          25     15  	16	 0.04    	0.04    	;     
	      26     7   	16 	 0.09    	0.12    	];
      
%    Loop data
Loop        = [11 12 15	19	18	16 0 0 0 0 0 0 0; 16 17	21	24	22 0 0 0 0 0 0 0 0; 13 14 26 25 23 24 21 17 18 19 15 12 20];
%Branch & possible parents
Parent =   [1   0   0   0   0   0;
            0   0   0   0   0   0;
            0   0   0   0   0   0;
            0   0   0   0   0   0;
            0   0   0   0   0   0;
            0   0   0   0   0   0;
            0   0   0   0   0   0;
            0   0   0   0   0   0;
            0   0   0   0   0   0;
            0   0   0   0   0   0;
            11  1   12  16  13  22;
            12  11  13  15  0   0 ;
            13  11  12  14  0   0;
            14  13  26  0   0   0;
            15  12  19  0   0   0;
            16  1   11  17  18  22;
            17  16  18  21  0   0;
            18  16  17  19  20  0;
            19  15  18  20  0   0;
            20  18  19  0   0   0;
            21  17  24  0   0   0;
            22  1   11  16  23  24;
            23  22  24  25  0   0;
            24  21  22  23  0   0;
            25  23  26   0  0   0;
            26  14  25  0   0   0];
root_branch = 1;                 %root branch

% Capacitor bank data
cap_bus=  [4 5 11 13 16 ];   %buses with capacitor banks
cap_list=[300 600 900 1200 1500 1800];    % capacitor bank values
%cap_list=[900 1200 1500 1800];    % capacitor bank values
% Initializations
Bus1(:,2:3)   = Bus1(:,2:3)/100;
Bus=Bus1;
cap_list=cap_list./100000;

sizeof_bus=size(Bus);
sizeof_line = size(Line);
sizeof_loop = size(Loop);
sizeof_parent = size(Parent);
sizeof_cap_bus=size(cap_bus);
sizeof_cap_list=size(cap_list);

Tie_list=[12 15 19; 18 17 21; 14 26 25];
sizeof_tie_list=size(Tie_list);
c = max(Tie_list);
d = sort(c);
phero_matrix_size = d(end);   %matrix size with maximum switch number
n_bus=phero_matrix_size;
cap_node=sizeof_cap_bus(2);

alpha =1; %5   %5          %1;	% Parameter representing the importance of trail
beta = 8 ;%10  %8         %5;	% Parameter representing the importance of visibility
rho = .5; %1   %1        %0.5;	% Evaporation
Q = 10 ; %10       %10;	% A constant

n_stage = sizeof_loop(1);    %Number of loops
ant_rec = n_stage;     % minimum no. of ants should be the number of tie options in stage1
ant_cap = sizeof_cap_bus(2);
Iter_max = 20;	% Maximum number of iterations

best_rec = zeros(Iter_max,n_stage);  %best path of each iteration
pheromone_rec = ones(phero_matrix_size,phero_matrix_size);    %pheromone content of paths from one stage to next
eta_rec = zeros(phero_matrix_size,phero_matrix_size);   %inverse of power loss corresponding to the above

ties=26; %total tie options to define pheromone and eta
pheromone_cap = ones(sizeof_cap_list(2)*cap_node,sizeof_cap_list(2)*cap_node,ties);    %pheromone content of paths from one stage to next
eta_cap = zeros(sizeof_cap_list(2)*cap_node,sizeof_cap_list(2)*cap_node,ties);   %inverse of power loss corresponding to the above

%***************************ITERATIONS*************************************
for iter_iter = 0:Iter_max-1      %iterations
    ITERATION=iter_iter
    t = cputime;
    rand_order1 = [];
    rand_position = [];
    Tie_list_stage1 = Tie_list(1,:);
    zero = find(Tie_list_stage1 == 0);
    for x = 1:ant_rec
        n_switch(x,1) = length(Tie_list_stage1)-size(zero,2); %size(zero,2) gives the total 0s in the tie_list_stage
    end    
    for i2 = 1:ceil(ant_rec/n_switch(1,1))
        rand_order1 = [rand_order1,randperm(n_switch(1,1))];
    end
    rand_order = rand_order1(1:ant_rec);
    for i2 = 1:length(rand_order)
        rand_position(i2) = Tie_list_stage1(rand_order(i2));
    end
    tabu_rec(:,1) = (rand_position(1:ant_rec))';    % the ants are placed randomly in each switch  
     best_cap_stage=zeros(n_stage,cap_node,ant_rec);
     
    %************************STAGES****************************************
    for iter_stage = 2:n_stage
        STAGE=iter_stage    
        Tie_list_stage  = Tie_list(iter_stage,:);
        sizeof_tie_list_stage=size(Tie_list_stage);
        Tie_list_stage1 = [];
        for x = 1:sizeof_tie_list_stage(2)
            if Tie_list_stage(x) ~= 0
                Tie_list_stage1=[Tie_list_stage1,Tie_list_stage(x)];
            end
        end
        sizeof_tie_list_stage1=size(Tie_list_stage1);       
        %****************************ANTS_RECONFIGURATION******************
        for iter_ant_rec = 1:ant_rec
        %iter_ant_rec
            switch_opened = tabu_rec(iter_ant_rec,1:(iter_stage-1)); %upto previous loop
            Tie_list_stage = [];
            Tie_list_stage2 = [];
            for x = 1:sizeof_tie_list_stage1(2)
                if Tie_list_stage1(x) ~= switch_opened
                    Tie_list_stage2 = [Tie_list_stage2,Tie_list_stage1(x)]; %avoid all previous ties in the present tie list options
                end
            end
            sizeof_tie_list_stage2 = size(Tie_list_stage2);
            ant_iter = 0;      %actual tie list size of the present ant in this stage
            Tie_open = zeros(1:n_stage);
            for j2 = 1:iter_stage-1
                Tie_open(j2) = tabu_rec(iter_ant_rec,j2);
            end
            Tie_list_stage5 = [];
            for x = 1:sizeof_tie_list_stage2(2) %check whether any node is left without supply
                Tie_open(iter_stage) = Tie_list_stage2(x);                
                [order_loop,continuous] = order_voltage(Bus,Line,Loop,Parent,Tie_open,iter_stage);
                if continuous == 1
                    ant_iter = ant_iter+1;  %changes from ant to ant in any stage
                    Tie_list_stage5 = [Tie_list_stage5,Tie_list_stage2(x)];
                end
            end
            Tie_list_stage = [];
            Tie_list_stage = Tie_list_stage5;
            sizeof_tie_list_stage = size(Tie_list_stage);
            n_switch(iter_ant_rec,iter_stage) = ant_iter;
            probability_rec = zeros(1,n_switch(iter_ant_rec,iter_stage)); 
            best_cap_switch=zeros(sizeof_tie_list(2),cap_node);
            %*********************TIE SWITCH*******************************
            for iter_switch = 1:n_switch(iter_ant_rec,iter_stage)  %each switch at a particular loop
            %iter_switch    
                Tie_open(iter_stage) = Tie_list_stage(iter_switch);               
                [tabu_cap,best_cap,Bus_cap,new_pheromone] = capacitor_place(Bus,Line,Loop,Parent,Tie_open,iter_stage,pheromone_cap,eta_cap);
                best_cap_switch(iter_switch,:)=best_cap;                
                [order_loop,continuous] = order_voltage(Bus_cap,Line,Loop,Parent,Tie_open,iter_stage);                
                [Total_powerloss_system,voltage_system] = powerloss_voltage(Bus_cap,Line,order_loop,iter_stage);
                eta_rec(switch_opened(end),Tie_list_stage(iter_switch)) = 1.0/Total_powerloss_system; %corresponding to previous loop-tie, current tie-option
                probability_rec(iter_switch) = (pheromone_rec(switch_opened(end),Tie_list_stage(iter_switch)))^alpha*(eta_rec(switch_opened(end),Tie_list_stage(iter_switch)))^beta;
            end
            probability_rec	= probability_rec/sum(probability_rec);
            pcum = cumsum(probability_rec);
            select = find(pcum >= rand);
            to_open = Tie_list_stage(select(1));
            tabu_rec(iter_ant_rec,iter_stage) = to_open; %tabu_list has information of each ants position at each stage corresponding to one iteration
            selected_tie=find(to_open==Tie_list(iter_stage,:));            
            best_cap_stage(iter_stage,:,iter_ant_rec)=best_cap_switch(selected_tie,:);            
        end  %iter_ant_rec
    end %iter_stage
    
    %calculating the minimum loss path and average loss in each iteration
    for i2 = 1:ant_rec    %total powerloss corresponding to each ant
        Tie_open = tabu_rec(i2,:);        
        [order_loop,continuous] = order_voltage(Bus_cap,Line,Loop,Parent,Tie_open,iter_stage);        
        [Total_powerloss_system,voltage_system] = powerloss_voltage(Bus_cap,Line,order_loop,iter_stage);
        total_loss(i2) = Total_powerloss_system;
    end
    loss_min(iter_iter+1) = min(total_loss);
    position = find(total_loss == loss_min(iter_iter+1));
    best(iter_iter+1,:) = tabu_rec(position(1),:); %best ant path of the iteration
    best_capacitor(iter_iter+1,:)=best_cap_stage(n_stage,:,position(1));
    loss_average(iter_iter+1) = mean(total_loss);
    delta_pheromone_rec	= zeros(phero_matrix_size,phero_matrix_size);
    ant1_loss(iter_iter+1) = total_loss(1);
    ant2_loss(iter_iter+1) = total_loss(2);
    ant3_loss(iter_iter+1) = total_loss(3);
   
    for i_1 = 1:ant_rec
        for j_1 = 1:(n_stage-1) %pheromone change 
            delta_pheromone_rec(tabu_rec(i_1,j_1),tabu_rec(i_1,j_1+1)) = delta_pheromone_rec(tabu_rec(i_1,j_1),tabu_rec(i_1,j_1+1))+Q/total_loss(i_1);
        end
    end
    ant = tabu_rec;    
    pheromone_rec = (1-rho).*pheromone_rec+delta_pheromone_rec;
    old_tabu_rec = tabu_rec;
    tabu_rec = zeros(ant_rec,n_stage);
    time= cputime-t;
    if best(iter_iter+1,:) == [19 17 26]
        count=count+1;
    end
end  %iter_iter
Total_loss_ant = total_loss
Solution_reconfiguration = best(Iter_max,:)
Solution_capacitor_placement=best_capacitor(Iter_max,:)
Solution_path_loss = loss_min(end)
Average_loss = loss_average(end) 
Iteration_time=time
best_reconfiguration=best
best_capacitor=best_capacitor
loss_min_iteration=loss_min'
count
%--------------------------------------------------------------------------
figure(1)
set(gcf,'Name','Ant Colony Optimization！！Figure of loos_min and loss_average','Color','w')
plot(loss_min,'r')
set(gca,'Color','w')
hold on
plot(loss_average,'k')
xlabel('Iterations')
ylabel('Min powerloss & Average powerloss')
title('Powerloss')

figure(2)    
set(gcf,'Name','Ant Colony Optimization！！Ant paths')
plot(ant(1,:),'y')
hold on
plot(ant(2,:),'m')
hold on
plot(ant(3,:),'c')
%hold on
%plot(ant(4,:),'r')
%hold on
%plot(ant(5,:),'g')
%hold on
%plot(ant(6,:),'b')
%hold on
%plot(ant(7,:),'k')
%hold on
xlabel('Stages')
ylabel('Tie Switches')
title('Ant paths')

Original=[0.9951;0.9951;0.9951;0.9912;0.9888;0.9861;0.9850;0.9777;0.9710;0.9770;0.9710;0.9690;0.9944;0.9950;0.9915;0.9910];
%Original=Bus(:,4)';
New=voltage_system(:,2)';
figure(3)
set(gcf,'Name','Voltage_Profile！！Original ang New Voltages','Color','w')
plot(Original,'k','LineWidth',1)
set(gca,'Color','w')
hold on
plot(New,'r','LineWidth',2)
axis([1 16 0.968 1])
xlabel('Nodes')
ylabel('Voltage')
title('Voltage_Profile')
