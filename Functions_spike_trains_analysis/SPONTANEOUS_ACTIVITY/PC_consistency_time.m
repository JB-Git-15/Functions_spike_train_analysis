function [Prods,Mean_s,STD_s,h] = PC_consistency_time(filename,tWin,N_jitt,minRate,plot_or_not)

%%% This function is intended to chop the data and verify the consistency
%%% of the overlap between the PC of the first part of the data and the
%%% second. The overlap = |PC_n(1st part).*PC_n(second part)|
%%% The first and the second part are defined cutting the data into parts
%%% and then permutting this parts and taking the first and the second part

%%% V1 in a first version, we cut the data into 6 parts


load(filename)
[SpkCountMat,SpkCountMatC,iLow,Spk_C_Mat_z]  = SpkCountMat_Centered_and_normalized(filename,tWin,N_jitt,minRate);

 
Spk_C_Mat_z                                  = Spk_C_Mat_z(1:floor( size(Spk_C_Mat_z,1)/(2*N_jitt))*2*N_jitt,:);  % divisible par N jitt et par 2
N_cells                                      = size(Spk_C_Mat_z,2);

Num_parts_J                                  = size(Spk_C_Mat_z,1)/N_jitt;
Num_perm                                     = 100;% min([ 100  Num_parts_J ]); 
Prods                                        = zeros(Num_perm,N_cells);                                            

        for    n = 1 : Num_perm 
           
           Perms =   randperm(Num_parts_J); 
           ind_f =   Perms(1: Num_parts_J/2);
           ind_s =   Perms(Num_parts_J/2 +1:end);
           Indices_first_half = []; 
           Indices_sec_half   = []; 
           
           for nn = 1 : length(ind_f)
              Indices_first_half = [Indices_first_half ; ( (ind_f(nn)- 1)*N_jitt + 1  : ind_f(nn)*N_jitt)' ];
              Indices_sec_half   = [Indices_sec_half   ; ( (ind_s(nn)- 1)*N_jitt + 1  : ind_s(nn)*N_jitt)' ];
           end
           
               First_half =   Spk_C_Mat_z(Indices_first_half,:) ; 
               Secnd_half =   Spk_C_Mat_z(Indices_sec_half  ,:) ; 
           
                   Time_length                      =   size(First_half,1);
                   C_f                              =   First_half'*First_half/Time_length;
                   [PC_f,eigs,exp]                  =   pcacov(C_f); 
                   C_s                              =   Secnd_half'*Secnd_half/Time_length;
                   [PC_s,eigs,exp]                  =   pcacov(C_s);  
           
             for jj = 1 : N_cells
                    Prods(n,jj)                     =   abs(PC_f(:,jj)'*PC_s(:,jj));  
             end      
        end                                 
                    STD_s                           =   std(Prods);   % standard deviation of the products                          
                    Mean_s                          =   mean(Prods);   
        
        if plot_or_not
            Col = rand([N_cells,3]);
         h = figure('Color','white');
           subplot(2,1,1)
             hold on 
             for jjj = 1 : N_cells
                  for n_perms = 1 : size(Perms,1)
                     plot( [jjj ,jjj ] ,[Prods(n_perms,jjj) Prods(n_perms,jjj)] ,'.','Color',Col( jjj,:))
                  end
             end
             hold off  
             xlabel('PC')
             ylabel('|PC_n(first half).PC_n(sec half)|')
             title([filename,', 6 parts, ','tWin: ',num2str(tWin),', N jitt: ',num2str(N_jitt) ] )
           subplot(2,1,2)
                    hold on
                   for u = 1 : N_cells   
                    plot([ u u ], [STD_s(u)  STD_s(u)],'.','Color',Col(u,:))
                   end 
                   hold off
                    xlabel('PC')
                    ylabel('std')
        else h = 0;
        end
             
        
end