function New_order = Reorder_8_shanks(A)

N = 8;


     D  = [0,10,4,3,1,.5,.2,.01;     %  distance function  (can decrease more or less streeply, starting from the diagonal)
        10,0,10,4,3,1,.5,.2;
        4,10,0,10,4,3,1,.5;
        3,4,10,0,10,4,3,1;
        1,3,4,10,0,10,4,3;          
        .5,1,3,4,10,0,10,4;
        .2,.5,1,3,4,10,0,10;
        .01,.2,.5,1,3,4,10,0];
 
 
All_perms               = perms([ 1 : N  ] );
Number_of_perms         = size(All_perms,1);
Value_of_cost_function  = zeros(Number_of_perms,1);


for u = 1 : Number_of_perms
    B_temp                        = A(All_perms(u,:), All_perms(u,:));    % one possible realisation 
    Value_of_cost_function(u,1)   = sum(sum(D.*B_temp));
end

    
[val,ind]       = max(Value_of_cost_function);

New_order       = All_perms(ind(1),:);    

 

end