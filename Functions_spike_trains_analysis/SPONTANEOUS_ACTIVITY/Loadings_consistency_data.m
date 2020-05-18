function [Loadings_all_perms,h] = Loadings_consistency_data(filename,N_PC,tWin,N_jitt,minRate,plot_or_not)

%  Mean_s,STD_s ? 

%%% This function is intended to chop the data and verify the consistency
%%% of the loadings  between the PC of the first part of the data and the
%%% second. 
%%% The first and the second part are defined cutting the data into parts
%%% and then permutting this parts and taking the first and the second part

load(filename)
[SpkCountMat,SpkCountMatC,iLow,Spk_C_Mat_z]  = SpkCountMat_Centered_and_normalized(filename,tWin,N_jitt,minRate);

 
Spk_C_Mat_z                                  = Spk_C_Mat_z(1:floor( size(Spk_C_Mat_z,1)/(2*N_jitt))*2*N_jitt,:);  % divisible par N jitt et par 2
N_cells                                      = size(Spk_C_Mat_z,2);

Num_parts_J                                  = size(Spk_C_Mat_z,1)/N_jitt;
Num_perm                                     = 100;% min([ 100  Num_parts_J ]); 
Loadings_all_perms                           = zeros(Num_perm,N_cells);                                            

%%%%  First determine the loadings with all the data
            Time_length                      =   size(Spk_C_Mat_z,1);
            C_a                              =   Spk_C_Mat_z'*Spk_C_Mat_z/Time_length;
            [PC_a,eigs,exp]                  =   pcacov(C_a); 
            loadings_a                       =   PC_a(:,N_PC);
            [val,ind_loadings]               =   sort(loadings_a);                         
            %%%% 

 
      for    n = 1 : Num_perm 
           
           Perms =   randperm(Num_parts_J); 
           ind_f =   Perms(1: Num_parts_J/2);
           %ind_s =   Perms(Num_parts_J/2 +1:end);
           Indices_first_half = []; 
           %Indices_sec_half   = []; 
           
           for nn = 1 : length(ind_f)
              Indices_first_half = [Indices_first_half ; ( (ind_f(nn)- 1)*N_jitt + 1  : ind_f(nn)*N_jitt)' ];
             % Indices_sec_half  = [Indices_sec_half   ; ( (ind_s(nn)- 1)*N_jitt + 1  : ind_s(nn)*N_jitt)' ];
           end
           
               First_half =   Spk_C_Mat_z(Indices_first_half,:) ; 
              % Secnd_half =   Spk_C_Mat_z(Indices_sec_half  ,:) ; 
           
                   Time_length                      =   size(First_half,1);
                   C_f                              =   First_half'*First_half/Time_length;
                   [PC_f,eigs,exp]                  =   pcacov(C_f); 
                %                    C_s                              =   Secnd_half'*Secnd_half/Time_length;
                %                    [PC_s,eigs,exp]                  =   pcacov(C_s);  
           
             
                       Loadings_all_perms(n,:)      =    PC_f(ind_loadings,N_PC)'*sign(PC_f(:,N_PC)'*loadings_a) ;
                 
        end                                 
                  %  STD_s                           =   std(Prods);   % standard deviation of the products                          
                  %  Mean_s                          =   mean(Prods);   
        
        if plot_or_not
            %Col = rand([N_cells,3]);
         h = figure('Color','white');
           
             hold on
                     plot(val,'r')
             for jjj = 1 : N_cells
                  for n_perms = 1 : size(Loadings_all_perms,1)
                     plot( [jjj ,jjj ] ,[Loadings_all_perms(n_perms,jjj) Loadings_all_perms(n_perms,jjj)] ,'g.' )
                  end
                     plot([jjj jjj],  2*[ std(Loadings_all_perms(:,jjj)) std(Loadings_all_perms(:,jjj))  ],'.m' ) 
                     plot([jjj jjj], -2*[ std(Loadings_all_perms(:,jjj)) std(Loadings_all_perms(:,jjj))  ],'.m' ) 

             end
                  plot(val,'r')
             hold off  
             xlabel('PC')
             ylabel('loadings consistency')
             title([filename(1 : end - 4),', PC : ',num2str(N_PC),', Perm :',num2str(Num_perm),', tWin: ',num2str(tWin),', N jitt: ',num2str(N_jitt) ] )
          
        else h = 0;
        end
             
        
end