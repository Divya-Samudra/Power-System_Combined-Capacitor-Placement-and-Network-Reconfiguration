function [Total_powerloss_system,voltage_system] = powerloss_voltage(Bus,Line,order_loop,iter_stage);
sizeof_order_loop = size(order_loop);
sizeof_line = size(Line);
sizeof_bus = size(Bus);
sizeof_order_loop = size(order_loop);

Pi_order_loop=[];
Qi_order_loop=[];
Vi_order_loop=[];
Ploss_order_loop = [];
Totalloss_order_loop = [];
Ploss_order_loop=zeros(sizeof_order_loop(1),sizeof_order_loop(2));
Totalloss_order_loop=zeros(sizeof_order_loop(1),sizeof_order_loop(2));
Pi_order_loop=zeros(sizeof_order_loop(1),sizeof_order_loop(2));
Qi_order_loop=zeros(sizeof_order_loop(1),sizeof_order_loop(2));
Vi_order_loop=zeros(sizeof_order_loop(1),sizeof_order_loop(2));
R=zeros(sizeof_order_loop(1),sizeof_order_loop(2));
X=zeros(sizeof_order_loop(1),sizeof_order_loop(2));
%---------------bus order
for x = 1:iter_stage
    for i1=1:sizeof_order_loop(1)
        zero=find(order_loop(i1,:)==0);
        if zero>0
            order_loop_now=order_loop(i1,1:zero(1)-1);
        else
            order_loop_now=order_loop(i1,:);
        end       
        sizeof_order_loop_now=size(order_loop_now);    
        for j1=2:sizeof_order_loop_now(2)-1        
            x1=find(order_loop(i1,sizeof_order_loop_now(2)-j1)==Line(:,1));
            x2=find(order_loop(i1,sizeof_order_loop_now(2)-j1+1)==Line(:,1));
            x3=find(Line(x1,2)==Line(x2,2:3));
            if x3>0
                a=Line(x1,2);            
            else
                a=Line(x1,3);            
            end
            order_loop_node(i1,sizeof_order_loop_now(2)-j1+1)=a;        
            if j1==2
                x1=find(order_loop(i1,1)==Line(:,1));
                x2=find(order_loop(i1,2)==Line(:,1));
                x3=find(Line(x1,2)==Line(x2,2:3));
                if x3>0
                    b=Line(x1,3);                
                else
                    b=Line(x1,2);               
                end
                order_loop_node(i1,1)=b;            
            end
        end
        order_loop_node(i1,sizeof_order_loop_now(2))=1; 
        if sizeof_order_loop_now(2)==2
            x1=find(order_loop(i1,1)==Line(:,1));
                x2=find(order_loop(i1,2)==Line(:,1));
                x3=find(Line(x1,2)==Line(x2,2:3));
                if x3>0
                    b=Line(x1,3);                
                else
                    b=Line(x1,2);               
                end
                order_loop_node(i1,1)=b;   
        end
    end
    sizeof_order_loop_node=size(order_loop_node);
    
    Pl=zeros(sizeof_order_loop_node(1),sizeof_order_loop_node(2));
    Ql=zeros(sizeof_order_loop_node(1),sizeof_order_loop_node(2));
    %-----------------------------------------------
         % impedances
      for i1=1:sizeof_order_loop(1)
        for j1 = 1:sizeof_order_loop(2)
            row_order_loop = order_loop(i1,:);
            pos_order_loop = find(row_order_loop == 0);
        end
        if pos_order_loop > 0% in case of last 0s row elements
            z1(i1) = pos_order_loop(1)-1;
        else
            z1(i1) = sizeof_order_loop(2);
        end        
        for k1 = 1:z1(i1) % actual branch size of each path is z1  
                %R, X from line input
                % Pl, Ql for every node from Bus input
                h=find(Line(:,1) == order_loop(i1,k1));
                    R(i1,k1) = Line(h,4);
                    X(i1,k1) = Line(h,5);
        end
      end
      %load real and reactive powers
      for i1=1:sizeof_order_loop_node(1)
        for j1 = 1:sizeof_order_loop_node(2)
            row_order_loop_node = order_loop_node(i1,:);
            pos_order_loop_node = find(row_order_loop_node == 0);
        end
        if pos_order_loop_node > 0% in case of last 0s row elements
            z2(i1) = pos_order_loop_node(1)-1;
        else
            z2(i1) = sizeof_order_loop_node(2);
        end        
        for k1 = 1:z2(i1) % actual branch size of each path is z1               
                % Pl, Ql for every node from Bus input
                h=find(Line(:,1) == order_loop_node(i1,k1));                                 
            Pl(i1,k1)=Bus(order_loop_node(i1,k1),2);
            Ql(i1,k1) = Bus(order_loop_node(i1,k1),3);
        end
      end
      % dist flow equations
      for i1=1:sizeof_order_loop_node(1)
        for j1 = 1:sizeof_order_loop_node(2)
            row_order_loop_node = order_loop_node(i1,:);
            pos_order_loop_node = find(row_order_loop_node == 0);
        end
        if pos_order_loop_node > 0% in case of last 0s row elements
            z3(i1) = pos_order_loop_node(1)-1;
        else
            z3(i1) = sizeof_order_loop_node(2);
        end 
        Vi_order_loop(i1,1)=Bus(order_loop_node(i1,1),4);
        Pi_order_loop(i1,1)=0;
        Qi_order_loop(i1,1)=0;
        for k1 = 2:z3(i1) 
            Pi_order_loop(i1,k1) = Pi_order_loop(i1,k1-1)+(R(i1,k1-1)*(((Pi_order_loop(i1,k1-1)+Pl(i1,k1-1))^2)+((Qi_order_loop(i1,k1-1)+Ql(i1,k1-1))^2))/Vi_order_loop(i1,k1-1)^2)+Pl(i1,k1-1);
            Qi_order_loop(i1,k1) = Qi_order_loop(i1,k1-1)+(X(i1,k1-1)*(((Pi_order_loop(i1,k1-1)+Pl(i1,k1-1))^2)+((Qi_order_loop(i1,k1-1)+Ql(i1,k1-1))^2))/Vi_order_loop(i1,k1-1)^2)+Ql(i1,k1-1);
            Vi_order_loop(i1,k1) = ((Vi_order_loop(i1,k1-1)^2)+(2*(R(i1,k1-1)*(Pi_order_loop(i1,k1-1)+Pl(i1,k1-1))+X(i1,k1-1)*(Qi_order_loop(i1,k1-1)+Ql(i1,k1-1))))+(((R(i1,k1-1)^2)+(X(i1,k1-1)^2))*((((Pi_order_loop(i1,k1-1)+Pl(i1,k1-1))^2)+((Qi_order_loop(i1,k1-1)+Ql(i1,k1-1))^2))/(Vi_order_loop(i1,k1-1)^2))))^(0.5);               
        end
      end
      loss_order_loop = 0;
    for i1 = 1:sizeof_order_loop(1)
        Totalloss_order_loop(i1) = 0.0;
        for j1 = 2:z1(i1)
            Ploss_order_loop(i1,j1) = R(i1,j1-1)*((Pi_order_loop(i1,j1)^2)+(Qi_order_loop(i1,j1)^2))/Vi_order_loop(i1,j1)^2;
            Totalloss_order_loop(i1) = Totalloss_order_loop(i1)+Ploss_order_loop(i1,j1);
        end
        loss_order_loop = loss_order_loop + Totalloss_order_loop(i1);
    end     
   %---------------------------------------------------------------------- 
end
%Total_powerloss_system = loss_order_loop.*100000;
%-----------------------------------------------
% voltage profile
count_volt=zeros(sizeof_bus(1),1);
voltage_profile=zeros(sizeof_bus(1),2);
voltage_profile(:,1)=Bus(:,1);
voltage_system=zeros(sizeof_bus(1),2);
voltage_system(:,1)=Bus(:,1);
for i1=1:sizeof_order_loop_node(1)
    for j1=1:sizeof_order_loop_node(2)
        voltage_value= find(voltage_profile(:,1)==order_loop_node(i1,j1));
        if voltage_value>0
            count_volt(voltage_value)=count_volt(voltage_value)+1;
            voltage_profile(voltage_value,2)=voltage_profile(voltage_value,2)+(Vi_order_loop(i1,j1)^2);
        end
    end
end
for i1=1:sizeof_bus(1)
    if count_volt(i1)~=0
    voltage_system(i1,2)=(voltage_profile(i1,2)/count_volt(i1))^(0.5);
    end
end
voltage_system(2,2)=voltage_system(1,2);
voltage_system(3,2)=voltage_system(1,2);

%---------------------------------------------------
%order_loop
%order_loop_node
%R
%X
%Pl
%Ql
%Pi_order_loop    
%Qi_order_loop        
%Vi_order_loop 
%Ploss_order_loop
Total_powerloss_system = loss_order_loop.*100000;  
%voltage_system
        
        
        
 