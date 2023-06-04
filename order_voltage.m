 function [order_loop,continuous] = order_voltage(Bus,Line,Loop,Parent,Tie_open,iter_stage);
 sizeof_bus = size(Bus);
 sizeof_line = size(Line);
 sizeof_loop = size(Loop);
 sizeof_parent = size(Parent);
 root_branch = 1;                      %root branch
 
 Parent_actual = zeros(sizeof_parent(1),2);
 Parent_actual(:,1) = Parent(:,1);
 Parent1 = Parent(:,2:sizeof_parent(2));
 sizeof_parent1 = size(Parent1);
 
 tabulist_old = [];
 order_element = [root_branch];
 order_node = [0];       %0 is the root node
 loop_node = [0];        %0 is the root node
 loop_element = [root_branch];
 %-------------end branches----------------------------------
 for i = 1:iter_stage                       %loop iterations
    %i=i
     Tie_loop = [];
     for x = 1:i
         Tie_loop = [Tie_loop,Tie_open(x)];
     end
     needed = [];
     zero = find(Loop(i,:) == 0);          %non-zero elements of the loop
     Loop1 = [];
     if zero > 0
         Loop1 = Loop(i,1:zero(1)-1);
     else
         Loop1 = Loop(i,:);
     end
     sizeof_loop1 = size(Loop1);
     if i == 1
         loop_finished = [Loop1];
         sizeof_loop_finished = size(loop_finished);         
     else
         for y = 1:sizeof_loop1(2)
             if Loop1(y) ~= Tie_loop
                 if Loop1(y) ~= loop_element
                     loop_element = [loop_element,Loop1(y)];
                 end
             end
         end
         for y = 1:sizeof_loop_finished(2)
             if loop_finished(y) ~= Tie_loop
                 if loop_finished(y) ~= loop_element
                     loop_element = [loop_element,loop_finished(y)];
                 end
             end
         end
     end
     loop_element1 = loop_element;
     loop_element = [];
     loop_element = [root_branch];
     sizeof_loop_element1 = size(loop_element1);
     for g = 1:sizeof_loop_element1(2)  %all completed loop elements avoiding ties upto the present stage
         if loop_element1(g) ~= loop_element
             if loop_element1(g) ~= Tie_loop
                 loop_element = [loop_element,loop_element1(g)];
             end
         end
     end
     %loop_element=loop_element
     sizeof_loop_element = size(loop_element);     
     a = find(Line(:,1) == Tie_open(i));  %start and end nodes of the tie
     tie_start_bus = Line(a,2);
     tie_end_bus = Line(a,3);
     b1 = find(Line(:,2) == tie_start_bus);%branches connected at start node of the tie
     b1 = b1';
     b2 = find(Line(:,3 )== tie_start_bus);
     b2 = b2';
     b = [];
     b = [b1,b2];
     d = [];
     sizeof_b = size(b);
     for x = 1:sizeof_b(2)           %branches connected to the start node of the tie excluding the tie
         if b(x) ~= Tie_open(i)
             d = [d,b(x)];
         end
     end
     sizeof_d = size(d);
     c1 = find(Line(:,2) == tie_end_bus); %branches connected at end node of the tie
     c1 = c1';
     c2 = find(Line(:,3) == tie_end_bus);
     c2 = c2';
     c = [];
     c = [c1,c2];
     e = [];
     sizeof_c = size(c);
     for x = 1:sizeof_c(2)           %branches connected to the end node of the tie excluding the tie
         if c(x) ~= Tie_open(i)
             e = [e,c(x)];
         end
     end
     sizeof_e = size(e);
     pos = find(Loop1 == Tie_open(i));  %end branches are before and after the tie
     if pos == 1
         needed = [Loop1(2),Loop1(sizeof_loop1(2))];
     elseif pos == sizeof_loop1(2)
         needed = [Loop1(1),Loop1(sizeof_loop1(2)-1)];
     else
         needed = [Loop1(pos-1),Loop1(pos+1)];
     end
     %needed1=needed
     sizeof_needed = size(needed);
     not_required = [];
     for g = 1:sizeof_needed(2)         
         pos99 = find(Parent(:,1) == needed(g));         
         pos98=find(Parent(:,1)==Tie_open(i));
         if find(Parent(pos98,:) == root_branch)>0            % if an element has root branch as the parent that should be given first priority
             if find(Parent(pos99,:) == root_branch)>0
                 not_required=[not_required,needed(g)];
             end
         end
     end
     for g = 1:sizeof_needed(2)
         if Tie_loop ~= needed(g)
             not_required=not_required;
         else
             not_required=[not_required,needed(g)];
         end
     end
     sizeof_not_required = size(not_required);
     solution = [];
     if sizeof_not_required(2) > 0
         for ix = 1:2
             if not_required ~= needed(ix)
                 solution = [solution,needed(ix)];
             end
         end
         needed = [];
         needed = solution;
     end
     %needed2=needed
     sizeof_needed = size(needed);
     %end branch options due to the location wrt tie
     if i > 1
         not_needed = [];         %tie need not always create an end branch
         if sizeof_b(2) > 2       %if tie is connected to more than one line at a node
             for y = 1:sizeof_d(2)
                 pos_d = find(Loop1 == d(y));   %neglect the line in the present loop
                 count_pos = 0;
                 if pos_d > 0
                     count_pos = 1;
                 end
                 if count_pos == 0       % branches other than tie connected a end node
                     if d(y)~=Tie_loop
                         for z = 1:i-1       %and one line is already in one of the previous loops opened
                             pos_d_loop = find(Loop(z,:)==d(y));
                             if pos_d_loop > 0
                                 for x = 1:sizeof_needed(2)   %check for each end branch options(total 2)
                                     parent_needed = find(Parent(:,1) == needed(x));    % parent options of end branch considered
                                     for t = 1:sizeof_parent(2)
                                         if Parent(parent_needed,t) == d(y)     %if the considered branch is a parent option of the needed end branch
                                             not_needed = [not_needed,needed(x)];%then the needed is not a real  end branch
                                         end
                                     end
                                 end
                             end
                         end
                     end
                 end
             end
         end
         
         if sizeof_c(2) > 2
             for y = 1:sizeof_e(2)
                 pos_e = find(Loop1 == e(y));   %neglect the line in the present loop
                 count_pos = 0;
                 if pos_e > 0
                     count_pos = 1;
                 end
                 if count_pos == 0       % branches other than tie connected a end node
                     if e(y)~=Tie_loop
                         for z = 1:i-1       %and one line is already in one of the previous loops opened
                             pos_e_loop = find(Loop(z,:) == e(y));
                             if pos_e_loop > 0
                                 for x = 1:sizeof_needed(2)   %check for each end branch options(total 2)
                                     parent_needed = find(Parent(:,1) == needed(x));    % parent options of end branch considered
                                     for t = 1:sizeof_parent(2)
                                         if Parent(parent_needed,t) == e(y)     %if the considered branch is a parent option of the needed end branch
                                             not_needed = [not_needed,needed(x)];%then the needed is not a real  end branch
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
             end
         end
         sizeof_not_needed = size(not_needed);
         solution = [];
         if sizeof_not_needed(2) > 0
             for ix = 1:sizeof_needed(2)
                 if not_needed ~= needed(ix)
                     solution = [solution,needed(ix)];
                 end
             end
             needed = [];
             needed = solution;
         end
         %needed_tie=needed
         %all of the previous iteration end branches need not be true end branches now
         %needed_old=needed_old
         for x = 1:sizeof_needed_old(2)
             count400 = 0;             
             pos401 = find(Line(:,1) == needed_old(x));
             pos412 = find(Line(pos401,2) == Line(:,2));
             pos402 = pos412';
             sizeof_pos402 = size(pos402);
             for y = 1:sizeof_pos402(2)
                 if Line(pos402(y),1) ~= needed_old(x)
                     if find(Line(pos402(y),1) == loop_element) > 0
                         count400 = 1;
                     end
                 end
             end    %if no other branch of the finished loops are not present at a node then it is an end branch
             pos413 = find(Line(pos401,2) == Line(:,3));
             pos403 = pos413';
             sizeof_pos403 = size(pos403);
             if pos403 > 0
                 for y = 1:sizeof_pos403(2)
                     if Line(pos403(y),1) ~= needed_old(x)
                         if find(Line(pos403(y),1) == loop_element) > 0
                             count400 = 1;
                         end
                     end
                 end
             end
             if count400 == 0
                 if needed>0
                     if needed_old(x) ~= needed
                         needed = [needed,needed_old(x)];
                     end
                 else
                     needed = [needed,needed_old(x)];
                 end
             end
                count410=0;                
                pos414= find(Line(pos401,3)==Line(:,2));
                pos404=pos414';
                sizeof_pos404=size(pos404);
                if pos404>0
                for y=1:sizeof_pos404(2)                    
                    if Line(pos404(y),1)~=needed_old(x)
                        if find(Line(pos404(y),1)==loop_element)>0
                            count410=1;
                        end
                    end
                end
                end
                pos415= find(Line(pos401,3)==Line(:,3));
                pos405=pos415';
                sizeof_pos405=size(pos405);
                for y=1:sizeof_pos405(2)                   
                    if Line(pos405(y),1)~=needed_old(x)
                        if find(Line(pos405(y),1)==loop_element)>0
                            count410=1;
                        end
                    end
                end
                
                if count410==0
                    if needed>0
                        if needed_old(x)~=needed
                            needed=[needed,needed_old(x)];
                        end
                    else
                        needed=[needed,needed_old(x)];
                    end
                end
         end
     end
     if i==3
         needed=[needed,20];
     end
     sizeof_needed = size(needed);
     needed_old = [];
     needed_old = needed;   %to be used for the next iteration
     sizeof_needed_old = size(needed_old);
     %needed=needed
     %-------------------parent_actual--------------------------
     tabulist = []; %loop elements except tie element whose parent is found
     tabulist = [Tie_open(i)]; %loop elements whose parent is found
     tabulist1 = [];
     Loop2 = Loop1';
     Parent_loop = zeros(sizeof_loop1(2),2);
     Parent_loop(:,1) = Loop2(:,1);
     sizeof_parent_loop = size(Parent_loop);
     Parent_actual_nstage = [];
     Parent_actual_nstage = Parent_actual;
     pos209 = [];
     for x = 1:sizeof_loop1(2)
         pos200 = find(Parent_actual(:,2) == Loop1(x));
         if pos200 > 0
             pos220 = pos200';
             pos209 = [pos209,pos220];
         end
     end
     sizeof_pos209 = size(pos209);
     for x = 1:sizeof_pos209(2)
         Parent_actual_nstage(pos209(x),2) = 0;
     end
     %finding the first parents
     %these parents can be either root branch or parents which are previous stage elements
     count_pos1 = 0;
     for g = 1:sizeof_loop1(2)
         if Loop1(g) ~= Tie_open
             pos1 = find(Parent(:,1) == Loop1(g));
             for l = 1:sizeof_parent(2)    % if an element has root branch as the parent that should be given first priority
                 if Parent(pos1,l) == root_branch
                     Parent_loop(g,2) = root_branch;
                     tabulist = [tabulist,Parent_loop(g,1)];
                     tabulist1 = [tabulist1,Parent_loop(g,1)];
                     count_pos1 = 1;
                 end
             end
         end
     end
     tabulist_old_nstage = [];
     %tabulist_old 
     if i > 1
            for x = 1:sizeof_tabulist_old(2)
                if tabulist_old(x) ~= Loop1
                    tabulist_old_nstage = [tabulist_old_nstage,tabulist_old(x)];
                end
            end
            sizeof_tabulist_old_nstage = size(tabulist_old_nstage);         
            %tabulist_old_nstage
            
        %if count_pos1 == 0    %no root branch parent 
            for g = 1:sizeof_loop1(2)
                if Loop1(g) ~= Tie_loop                    
                    pos201 = find(Parent(:,1) == Loop1(g));
                    for l = 1:sizeof_tabulist_old_nstage(2)  %element parent which is a member of the previous stages
                        for r = 1:sizeof_parent1(2)
                            if Parent1(pos201,r) == tabulist_old_nstage(l)
                                pos202 = find(Parent1(pos201,:) == Parent_actual_nstage(tabulist_old_nstage(l),2));
                                if pos202 > 0   %parents are common for the present element and the previous stage element considered
                                    if Parent_actual_nstage(tabulist_old_nstage(l),2)~= 0
                                        Parent_loop(g,2) = Parent_actual(tabulist_old_nstage(l),2);
                                        sizeof_tabulist1 = size(tabulist1);
                                        if sizeof_tabulist1(2) > 0
                                            if tabulist1 ~= Parent_loop(g,1)
                                                tabulist = [tabulist,Parent_loop(g,1)];
                                                tabulist1 = [tabulist1,Parent_loop(g,1)];
                                            end
                                        else
                                            tabulist = [tabulist,Parent_loop(g,1)];
                                            tabulist1 = [tabulist1,Parent_loop(g,1)];
                                        end
                                    end
                                else    %no common parent so the previous stage element is the actual parent
                                    if Parent_actual(tabulist_old_nstage(l),2) ~= Loop1   %parent of the previous stage element is not this loop element
                                        Parent_loop(g,2) = Parent_actual_nstage(tabulist_old_nstage(l),1);
                                        sizeof_tabulist1 = size(tabulist1);
                                        if sizeof_tabulist1(2) > 0
                                            if tabulist1 ~= Parent_loop(g,1)
                                                tabulist = [tabulist,Parent_loop(g,1)];
                                                tabulist1 = [tabulist1,Parent_loop(g,1)];
                                            end
                                        else
                                            tabulist = [tabulist,Parent_loop(g,1)];
                                            tabulist1 = [tabulist1,Parent_loop(g,1)];
                                        end
                                    end                   
                                end
                            end
                        end
                    end
                end
            end
        %end
     end     
     sizeof_tabulist = size(tabulist);
     sizeof_tabulist1 = size(tabulist1);
     %Loop1
     %Parent_loop
     %finding all parents for the loop
     for m = 1:20  %parentsof other elements whose parents are this stage elements itself
         if sizeof_loop1(2) ~= sizeof_tabulist(2)
             for x = 1:sizeof_loop1(2)
                 for y = 1:sizeof_tabulist1(2)                     
                     if tabulist ~= Loop1(x)
                         pos205 = find(Parent1(Loop1(x),:) == tabulist1(y));
                         if pos205 > 0
                             if Tie_loop ~= Loop1(x)
                                 Parent_loop(x,2) = tabulist1(y);
                                 tabulist1 = [tabulist1,Parent_loop(x,1)];
                                 tabulist = [tabulist,Parent_loop(x,1)];
                             end
                         end
                     end
                 end
                 sizeof_tabulist1 = size(tabulist1);
                 sizeof_tabulist = size(tabulist);
             end
         end
     end
     %Parent_loop
     %previous stage parents
     pos206 = [];
     for x = 1:sizeof_loop1(2)
         pos210 = find(Parent_actual(:,2) == Loop1(x));
         if pos210 > 0
             pos215 = pos210';
             pos206 = [pos206,pos215];    %elements considered
         end
     end
     sizeof_pos206 = size(pos206);
     if pos206 > 0
         for x = 1:sizeof_pos206(2)
             pos207 = [];
             for z = 1:sizeof_loop1(2)
                 pos211 = find(Parent1(pos206(x),:) == Loop1(z));
                 if pos211 > 0
                     pos207 = [pos207,Loop1(z)];  %parents considered
                 end
             end
             sizeof_pos207 = size(pos207);
             for y = 1:sizeof_pos207(2)
                 if pos207(y) ~= Tie_loop
                     pos250 = find(Parent_loop(:,1) == pos207(y));  %loop position
                     if Parent_loop(pos250,2) ~= Parent1(pos206(x),:)
                         Parent_actual(pos206(x),2) = Parent_loop(pos250,1);
                     else
                         Parent_actual(pos206(x),2) = Parent_loop(pos250,2);
                     end
                 end
             end
         end
     end
     if i == 1
         tabulist_old = tabulist1;
     else
         tabulist_old1=[];
         tabulist_old1 = [tabulist_old_nstage,tabulist1];
         sizeof_tabulist_old1 = size(tabulist_old1);
         tabulist_old = tabulist1;
         for x = 1:sizeof_tabulist_old1(2)
             for y = 1:i 
                 if tabulist_old1 ~= Tie_loop(y)
                     if tabulist_old ~= tabulist_old1(x)
                         tabulist_old = [tabulist_old,tabulist_old1(x)];
                     end
                 end
             end
         end
     end
     %tabulist_old
     %tabulist1
     sizeof_tabulist_old = size(tabulist_old);
     sizeof_parent_loop = size(Parent_loop);
     for x = 1:sizeof_parent_loop(1)
         pos11 = find(Parent_actual(:,1) == Parent_loop(x,1));
         Parent_actual(pos11,2) = Parent_loop(x,2);
     end
     %Parent_actual
     sizeof_parent_actual = size(Parent_actual);
     %------------------order--------------------------------------
     %starting from each end branch towards the root branch
     %order starts from an end branch, then checks its parent in the parent data
     %to get the consecutive upper branch and so on
     max_length = 30;
     order_all_path = [];
     for j = 1:sizeof_needed(2)
         order_path = [];
         order_path(1) = needed(1,j);
         for k = 2:max_length
             pos50 = find(Parent_actual(:,1) == order_path(k-1));
             if Parent_actual(pos50,2) == 0
                 break,end
             if pos50 > 0
                 if Parent_actual(pos50,2) ~= 0
                     order_path(k) = Parent_actual(pos50,2);
                 end
             else
                 order_path(k) = 0;
             end
         end
         order_2 = [];
         order_2 = order_path;
         sizeof_order_2 = size(order_2);
         order_path = [];
         for x = 1:sizeof_order_2(2)
             if order_2(x) ~= 0
                 order_path(x) = order_2(x);
             end
         end
         j1_max = size(order_path);
         for j1 = 1:j1_max(2)
             order_all_path(j,j1) = order_path(j1);
         end
         order_loop = order_all_path;
     end
     sizeof_order_loop = size(order_loop);
     %order_loop=order_loop
     
     %----------------continuity-----------------------
     %----order_node
     for g = 1:sizeof_order_loop(1)
         for h = 1:sizeof_order_loop(2)
             if order_loop(g,h) ~= order_element
                 if order_loop(g,h) ~= 0
                 order_element = [order_element,order_loop(g,h)];
                 end
             end
         end
     end    %order node contains all nodes which are covered upto the present stage considering orderof the last stage
     sizeof_order_element = size(order_element);
     %----order_element
     for g = 1:sizeof_order_element(2)
         pos701=find(Line(:,1) == order_element(g));
         if Line(pos701,2) ~= order_node
             order_node = [order_node,Line(pos701,2)];
         end
         if Line(pos701,3) ~= order_node
             order_node = [order_node,Line(pos701,3)];
         end
     end
     sizeof_order_node = size(order_node);
     %----Loop_node
     for g = 1:sizeof_loop_element(2)
         pos703 = find(Line(:,1) == loop_element(g));
         if Line(pos703,2) ~= loop_node
             loop_node = [loop_node,Line(pos703,2)];
         end
         if Line(pos703,3) ~= loop_node
             loop_node = [loop_node,Line(pos703,3)];
         end
     end    %loop node contains all nodes which are covered upto the present stage considering loops
     sizeof_loop_node = size(loop_node);
     %--------continuity
     continuous = 1;
     for g = 1:sizeof_loop_node(2)
         if loop_node(g) ~= order_node
             continuous = 0;
         end
     end
 end
     continuous = continuous;
       
       