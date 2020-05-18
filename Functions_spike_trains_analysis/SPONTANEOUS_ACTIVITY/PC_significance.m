function [ p_vals] = PC_significance(filename,tWin,N_jitt,minRate,N_perms,plot_or_not)

%%%    This function computes the consistency of the PC_x direction : 
%%%   We cut the data in two parts. We extract the n th PC of the first part.
%%%   Then we take the other half and we compute the absolute value of the scalar product of the nth PC (first half) with nth PC second half
%%%   Then we do shuffling of the second part : in windows of size
%%%   N_jitt*tWin, we permute the spike counts so as to destroy
%%%   correlations among cells. We then compute the first PC of the shuffled data 
%%%   and make the product with the PC1 of the first half of the data.
%%%   Example of usage:  [  p_vals] = PC1_significance(filename,tWin,N_jitt,minRate,N_perms,plot_or_not)

load(filename)
[SpkCountMat,SpkCountMatC,iLow,Spk_C_Mat_z]  = SpkCountMat_Centered_and_normalized(filename,tWin,N_jitt,minRate);
%Spk_C_Mat_z                                 = Spk_C_Mat_z(1:floor( size(Spk_C_Mat_z,1)/N_jitt)*N_jitt,:);
First_half_C_z                               = Spk_C_Mat_z(1:floor(size(Spk_C_Mat_z,1)/2),:);
First_half_C_z                               = First_half_C_z(1:floor(size(First_half_C_z,1)/N_jitt)*N_jitt,:); % time length is divisible by N jitter
Secnd_half_C_z                               = Spk_C_Mat_z(1+floor(size(Spk_C_Mat_z,1)/2):end,:);
Secnd_half_C_z                               = Secnd_half_C_z(1:floor(size(Secnd_half_C_z,1)/N_jitt)*N_jitt,:);


N_cells                                      = size(Spk_C_Mat_z,2); 
Scalar_product_first_second_half             = zeros(N_cells,1);  % ith coordinate : scalar prod PCi (first half) . PCi (second half)
Scalar_product_first_second_half_shuffled    = zeros(N_cells,N_perms);
p_vals                                       = zeros(N_cells,1);


%%%%  First determine the PC s of the first half of the data
            Time_length_half                 =   size(First_half_C_z,1);
            C_f                              =   First_half_C_z'*First_half_C_z/Time_length_half;    % correlation matrix
            [PC_f,eigs,exp]                  =   pcacov(C_f); 
             
%%%%  Second determine the PC s of the second half of the data
            Time_sec_half                    =   size(Secnd_half_C_z,1);
            C_s                              =   Secnd_half_C_z'*Secnd_half_C_z/Time_sec_half;
            [PC_s,eigs,exp]                  =   pcacov(C_s); 
                        
%%%% Compute the absolute value of the scalar product between the first and
%%%% the second half, for all the PCs

            for n_PC = 1 : N_cells
                 Scalar_product_first_second_half(n_PC) = abs(PC_f(:,n_PC)'*PC_s(:,n_PC));
            end
%%%%  Do the shuffling of the second half

  for n_perm = 1 : N_perms
                        Secnd_half_C_z_permuted            = zeros(size(Secnd_half_C_z));
                        nn                                 = N_jitt;
    while nn<= Time_sec_half
                       temp                                = Secnd_half_C_z(nn-N_jitt+1:nn,:);
         for n_cell = 1 : N_cells
            Secnd_half_C_z_permuted(nn-N_jitt+1:nn,n_cell) = temp(randperm(N_jitt),n_cell);  
         end
                                                 nn        = nn + N_jitt;
    end
                                     C_shuffled            =   Secnd_half_C_z_permuted'*Secnd_half_C_z_permuted/Time_sec_half;
                                    [PC_sh,eigs,exp]       =   pcacov(C_shuffled); 
            for n_PC = 1 : N_cells
    Scalar_product_first_second_half_shuffled(n_PC,n_perm) = abs(PC_f(:,n_PC)'*PC_sh(:,n_PC));
            end
  end
            
%%%% Calculate the p value

    for n_cell = 1 : N_cells
       p_vals(n_cell) = length(find(abs(Scalar_product_first_second_half(n_cell))< abs(Scalar_product_first_second_half_shuffled(n_cell,:))))/N_perms ;    
    end     
   
         if plot_or_not
             figure('Color','white')
             subplot(211)
             hold on
              for n_cell = 1 : N_cells
                 plot( n_cell*ones(N_perms,1), Scalar_product_first_second_half_shuffled(n_cell,:) ,'g.' )     
                 plot( [n_cell  n_cell],[ Scalar_product_first_second_half(n_cell) Scalar_product_first_second_half(n_cell) ],'r+' )     

               end 
             hold off
               title([filename, '# shuffles : ', num2str(N_perms) ,' ,t_c = ',num2str(tWin),' (s), N_J =',num2str(N_jitt)] )
               xlabel('PC #')
               ylabel('|PC first.PC second half |')
             subplot(212)
             hold on
             for nn = 1 : N_cells
                 if p_vals(nn) == 0
                  plot( [nn nn ],[ p_vals(nn)  p_vals(nn)],'r.')
                 else
                  plot( [nn nn ],[ p_vals(nn)  p_vals(nn)],'.')                
                 end
             end 
             hold off
              xlabel('PC #')
              ylabel('p_{val} ')
         end

  
        
end