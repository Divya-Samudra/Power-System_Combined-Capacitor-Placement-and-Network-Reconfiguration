function [tabu_cap,best_cap,Bus_cap,new_pheromone] = capacitor_place(Bus,Line,Loop,Parent,Tie_open,iter_stage,pheromone_cap,eta_cap);
        
% Capacitor bank data
cap_bus=  [4 5 11 13 16 ];   %buses with capacitor banks
cap_list=[300 600 900 1200 1500 1800];    % capacitor bank values
%cap_list=[900 1200 1500 1800];

root_branch = 1;                 %root branch

% Initializations
cap_list=cap_list./100000;
sizeof_bus=size(Bus);
sizeof_line = size(Line);
sizeof_loop = size(Loop);
sizeof_parent = size(Parent);
sizeof_cap_bus=size(cap_bus);
sizeof_cap_list=size(cap_list); 
cap_node=sizeof_cap_bus(2);

c = max(Bus(:,1));
d = sort(c);
n_bus = d(end);   % maximum bus number


ant_cap = sizeof_cap_bus(2);    % no. of ants 
n_stage = sizeof_loop(1);    %Number of loops
Iter_max_cap = 10;	% Maximum number of iterations
alpha =1; %5   %5          %1;	% Parameter representing the importance of trail
beta = 8 ;%10  %8         %5;	% Parameter representing the importance of visibility
rho = .5; %1   %1        %0.5;	% Evaporation
Q = 10 ; %10       %10;	% A constant

best_cap = zeros(Iter_max_cap,cap_node);  %best path of each iteration
%pheromone_cap = ones(sizeof_cap_list(2)*cap_node,sizeof_cap_list(2)*cap_node);    %pheromone content of paths from one stage to next
%eta_cap = zeros(sizeof_cap_list(2)*cap_node,sizeof_cap_list(2)*cap_node);   %inverse of power loss corresponding to the above
tie_present=Tie_open(iter_stage);

req_node_old=[];
%*********************Loops********************************
for iter_stage_cap=1:iter_stage    
%iter_stage_cap         
    % nodes to which capacitors are to be connected
    req_node=[];
    Loop2=Loop(iter_stage_cap,:);
    Loop1=[];
    for x=1:sizeof_loop(2)
        if Loop2(x)~=0
            Loop1=[Loop1,Loop2(x)];
        end
    end
    sizeof_loop1=size(Loop1);
    for x=1:sizeof_loop1(2)        %nodes having capacitor bank for each loop
        cap=find(Loop1(x)==Line(:,1));
        cap1=find(Line(cap,2)==cap_bus);
        if cap1>0
            pos1=find(Line(cap,2)==req_node);
            if pos1>0
                req_node=[req_node];
            else
                req_node=[req_node,Line(cap,2)];
            end
        end
        cap2=find(Line(cap,3)==cap_bus);
        if cap2>0
            pos2=find(Line(cap,3)==req_node);
            if pos2>0
                req_node=[req_node];
            else
                req_node=[req_node,Line(cap,3)];
            end
        end
    end
    sizeof_req_node=size(req_node);
    % positions of required nodes in capacitor bank
    position=[];
    for i1=1:sizeof_req_node(2)
        pos=find(req_node(i1)==cap_bus);
        position=[position,pos];
    end
    %all previous loop compensated nodes
    if iter_stage_cap>1
        for i=1:sizeof_req_node(2)
            if req_node(i)~=req_node_old
                req_node_old=[req_node_old,req_node(i)];
            end
        end
    else
        req_node_old=req_node;
    end
    req_node=[];
    req_node=req_node_old; 
    sizeof_req_node=size(req_node);
%req_node
end     %iter_stage
%*********************Iterations***************************
for iter_iter_cap=0:Iter_max_cap-1
%iter_iter_cap
    tabu_cap=zeros(ant_cap,sizeof_cap_bus(2));
    % at required node 1 ants are placed randomly at each capacitor value
    rand_order1=[];
    for i2=1:ceil(ant_cap/sizeof_cap_list(2))  %randomly positioning each ant at the first node of a loop
        rand_order1=[rand_order1,randperm(sizeof_cap_list(2))];
    end
    rand_order=rand_order1(1:ant_cap);
    for i2=1:length(rand_order)
        rand_position(i2)=cap_list(rand_order(i2));
    end
    tabu_cap(:,1)=(rand_position(1:ant_cap))';
    Bus2=Bus;
    %**********************Nodes*******************************
    for iter_node=2:sizeof_req_node(2)
        pos_node=find(req_node(iter_node)==cap_bus);
        %***********************Ants*******************************
        for iter_ant=1:ant_cap
%iter_ant 
            pos_node_old=find(req_node(iter_node-1)==cap_bus);
            cap_opened=tabu_cap(iter_ant,1:iter_node-1); %already opened capacitor values of previous nodes             
            probability_cap=zeros(1,sizeof_cap_list(2));
            pos_cap_old=find(cap_opened(end)==cap_list);
            eta_pos1=pos_node_old*sizeof_cap_list(2)-(sizeof_cap_list(2)-pos_cap_old);
            %***********************Capacitor switches*****************
            for iter_cap=1:sizeof_cap_list(2)
                cap_opened(iter_node)=cap_list(iter_cap); 
                for i=1:iter_node
                    node=find(req_node(i)==Bus(:,1));
                    Bus2(node,3)=Bus(node,3)-cap_opened(iter_node);
                end
                pos_cap=find(cap_list(iter_cap)==cap_list);
                [order_loop,continuous] = order_voltage(Bus2,Line,Loop,Parent,Tie_open,iter_stage_cap);
                [Total_powerloss_system,voltage_system] = powerloss_voltage(Bus2,Line,order_loop,iter_stage_cap);
                eta_pos2=pos_node*sizeof_cap_list(2)-(sizeof_cap_list(2)-pos_cap);
                eta_cap(eta_pos1,eta_pos2,tie_present)=1.0/Total_powerloss_system;
                probability_cap(iter_cap)=(pheromone_cap(eta_pos1,eta_pos2,tie_present))^alpha*(eta_cap(eta_pos1,eta_pos2,tie_present))^beta;
            end %iter_cap
            probability_cap=probability_cap/sum(probability_cap);% normalized probability of each capacitor in the bus capacitor bank
            pcum=cumsum(probability_cap);
            select=find(pcum>=rand);
            to_add=cap_list(select(1));
            tabu_cap(iter_ant,iter_node)=to_add;
        end %iter_ant
    end %iter_node    
    Bus2=Bus;
    
    for i2=1:ant_cap  %calculating powerloss with the selected capacitor values at every nodes
        Bus3=Bus;
        for j2=1:sizeof_req_node(2)
            node_1=find(req_node(j2)==Bus(:,1));
            pos_2=find(req_node(j2)==cap_bus);
            Bus3(node_1,3)=Bus(node_1,3)-tabu_cap(i2,pos_2);
        end
        [order_loop,continuous] = order_voltage(Bus3,Line,Loop,Parent,Tie_open,iter_stage_cap);
        [Total_powerloss_system,voltage_system] = powerloss_voltage(Bus3,Line,order_loop,iter_stage_cap);
        total_loss(i2)=Total_powerloss_system;
    end
    loss_min(iter_iter_cap+1)=min(total_loss);
    cap_placement(:,:)=tabu_cap(:,:);
    best_position=find(total_loss==loss_min(iter_iter_cap+1));
    best(iter_iter_cap+1,:)=cap_placement(best_position,:);
    loss_average(iter_iter_cap+1)=mean(total_loss);
    delta_pheromone_cap=zeros(sizeof_cap_list(2)*cap_node,sizeof_cap_list(2)*cap_node);
    for i_1=1:ant_cap
        for j_1=1:(sizeof_req_node(2)-1)
            x1=find(req_node(j_1)==cap_bus);
            x2=find(cap_placement(i_1,x1)==cap_list);
            x3=find(req_node(j_1+1)==cap_bus);
            x4=find(cap_placement(i_1,x3)==cap_list);
            delta_pheromone_cap(x1*sizeof_cap_list(2)-(sizeof_cap_list(2)-x2),x3*sizeof_cap_list(2)-(sizeof_cap_list(2)-x4))=delta_pheromone_cap(x1*sizeof_cap_list(2)-(sizeof_cap_list(2)-x2),x3*sizeof_cap_list(2)-(sizeof_cap_list(2)-x4))+Q/total_loss(i_1);
        end
    end    
    pheromone_cap(:,:,tie_present)=(1-rho).*pheromone_cap(:,:,tie_present)+delta_pheromone_cap;
    tabu_cap=zeros(ant_cap,cap_node);
    %[Total_powerloss_system,voltage_system] = powerloss_voltage(Bus3,Line,order_loop,iter_stage);    
end %iter_iter
Bus4=Bus;
for j=1:sizeof_cap_bus(2)
    node_final=find(cap_bus(j)==Bus(:,1));                
    Bus4(node_final,3)=Bus(node_final,3)-best(Iter_max_cap,j);
end
[order_loop,continuous] = order_voltage(Bus4,Line,Loop,Parent,Tie_open,iter_stage);
[Total_powerloss_system,voltage_system] = powerloss_voltage(Bus4,Line,order_loop,iter_stage);
best_cap=best(Iter_max_cap,:);
Bus_cap=Bus4;
new_pheromone=pheromone_cap;

   