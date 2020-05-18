clear all
close all
clc



 N = 9; % %  dimension 
 
 A = [0,7,5,3,2,2,3,5,7;
      7,0,7,5,3,2,2,3,5;
      5,7,0,7,5,3,2,2,3;
      3,5,7,0,7,5,3,2,2;
      2,3,5,7,0,7,5,3,2;
      2,2,3,5,7,0,7,5,3;
      3,2,2,3,5,7,0,7,5;
      5,3,2,2,3,5,7,0,7;
      7,5,3,2,2,3,5,7,0];
  
  figure;
  imagesc(A);title('A')
  
P = randperm(N);
B = A(P,P); 
figure;
imagesc(B); title('B : A scrambled')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D =  [0,7,5,3,2,2,3,5,7;        % distance function  (can decrease more or less streeply, starting from the diagonal)
      7,0,7,5,3,2,2,3,5;        %  play with this function eventually and automatise it ...
      5,7,0,7,5,3,2,2,3;
      3,5,7,0,7,5,3,2,2;
      2,3,5,7,0,7,5,3,2;
      2,2,3,5,7,0,7,5,3;
      3,2,2,3,5,7,0,7,5;
      5,3,2,2,3,5,7,0,7;
      7,5,3,2,2,3,5,7,0];


All_perms               = perms([ 1 : N  ] );
Number_of_perms         = size(All_perms,1);
Value_of_cost_function  = zeros(Number_of_perms,1);


for u = 1 : Number_of_perms
    B_temp                        = B(All_perms(u,:), All_perms(u,:));    % one possible realisation 
    Value_of_cost_function(u,1)   = sum(sum(D.*B_temp));
end

    
[val,ind]       = max(Value_of_cost_function);

Retained_permut = All_perms(ind(1),:); 

B_new = B(Retained_permut,Retained_permut);

figure;
imagesc(B_new); title('A rescued')















