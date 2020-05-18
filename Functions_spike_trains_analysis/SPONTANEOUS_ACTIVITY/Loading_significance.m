function [corr_score_cell,loadings,p_vals,h_1,Pente] = Loading_significance(filename,N_PC,tWin,N_jitt,minRate,N_perms,plot_or_not)

%%%    This function computes the significance of the loadings : 
%%%   for each cell i, lets consider the PCA of the correlation matrix of
%%%   the whole dataset without this cell i .Then, we compute the correlation
%%%   between the N_PC score and the counts of the cell i. 
%%%   In addition, we calculate the pvalue associated with each cell X, as
%%%   the proportion of cells in the data (we exclude cell X from the data )
%%%   that have a loading as extreme as the one from cell X.
%%%   h_1 is the handle of the image
%%%   Example of usage:   [corr_score_cell,loadings,p_vals,h_1] = Loading_significance(filename,.7,5,.5)

load(filename)
[SpkCountMat,SpkCountMatC,iLow,Spk_C_Mat_z]  = SpkCountMat_Centered_and_normalized(filename,tWin,N_jitt,minRate);
Spk_C_Mat_z                                  = Spk_C_Mat_z(1:floor( size(Spk_C_Mat_z,1)/N_jitt)*N_jitt,:);
N_cells                                      = size(Spk_C_Mat_z,2); 
corr_score_cell                              = zeros(N_cells,1);  
corr_score_cell_perm                         = zeros(N_cells,N_perms);
p_vals                                       = zeros(N_cells,1);

%%%%  First determine the loadings with all the data
            Time_length                      =   size(Spk_C_Mat_z,1);
            C_a                              =   Spk_C_Mat_z'*Spk_C_Mat_z/Time_length;
            [PC_a,eigs,exp]                  =   pcacov(C_a); 
            loadings                         =   PC_a(:,N_PC);
             
%%%% Substract a cell at a time and calculate the score and the p values
        for n_cell = 1 : N_cells
            
            SpkCount_z_temp                  = Spk_C_Mat_z;
            Count_cell_z                     = SpkCount_z_temp(:,n_cell);
            SpkCount_z_temp(:,n_cell)        = [];       
            C                                = SpkCount_z_temp'*SpkCount_z_temp/Time_length; % corr matrix
            [PC,eigs_t,exp]                  = pcacov(C);   
            
            PC_a_temp                        = loadings;
            PC_a_temp(n_cell)                = [];
            PC_N                             = PC(:,N_PC);
            PC_N                             = PC_N*sign(PC_a_temp'*PC_N);     % the direction of PC1 is arbitrary
            Score                            = SpkCount_z_temp*PC_N;
            %%% Normalize the variance of the score, before computing the
            %%% correlation
            Score                            = Score/sqrt(eigs_t(N_PC));
            %%% Compute the correlation
              corr_score_cell(n_cell)        = sum(Score.*Count_cell_z)/Time_length;
                                       %  p_vals(n_cell)                 = length(find(  abs(
                                       %  loadings(n_cell))  <  abs(PC(:,1)) ))/ (N_cells - 1) ;                         %  old story 
    %%%%  Do the shuffling to assess the statistical significance of the correlation between the cell and the score.            
     
  for n_perm = 1 : N_perms
         Count_cell_z_perm                  = zeros(Time_length,1);
         nn                                 = N_jitt;
    while nn<= Time_length
        temp                                = Count_cell_z(nn-N_jitt+1:nn);
        Count_cell_z_perm(nn-N_jitt+1:nn)   = temp(randperm(N_jitt));  
        nn                                  = nn + N_jitt;
    end
    corr_score_cell_perm(n_cell,n_perm)     = sum(Score.*Count_cell_z_perm)/Time_length;
 end           
     p_vals(n_cell)                         = length(find(abs(corr_score_cell_perm(n_cell,:)) > abs(corr_score_cell(n_cell))))/N_perms;   

        end
        
%%%%% Do the fitting corr = f(loading), save the slope        
     
        
                                               Y      = corr_score_cell;   
                                               X      = loadings; 
        
                                               Pente  = X\Y;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
         
        [PC_a_sorted, ind] = sort(loadings);
        if plot_or_not
           h_1 = figure('Color','white');
           
            % subplot(211)
            
         a1   = axes('position',[0.3 0.3    0.4    .4]);
             hold on
              for n_cell = 1 : N_cells
                  
                   if  ( corr_score_cell(ind(n_cell)) < max( corr_score_cell_perm(ind(n_cell),:)) ) &&   ( corr_score_cell(ind(n_cell)) > min( corr_score_cell_perm(ind(n_cell),:)) )
                 
              %   plot(repmat(PC_a_sorted(n_cell),[N_perms  1 ]),         corr_score_cell_perm(ind(n_cell),:),'g.' )      
                 plot(a1,[      PC_a_sorted(n_cell) PC_a_sorted(n_cell) ], [corr_score_cell(ind(n_cell))   corr_score_cell(ind(n_cell)) ],'o','Color',[0 0 0 ],'markerfacecolor', [1 1 1]) 
                  
                   else
                 plot(a1,[      PC_a_sorted(n_cell) PC_a_sorted(n_cell) ], [corr_score_cell(ind(n_cell))   corr_score_cell(ind(n_cell)) ],'o','Color',[0 0 0 ],'markerfacecolor', [0 0 0 ]) 
  
                   end
              end 
             hold off
               title(['  PC #', num2str(N_PC)] ,  'FontSize',12, 'fontWeight','Demi'  )  % filename, '# shuffles : ', num2str(N_perms) ,
               xlabel('loadings','FontSize',12, 'fontWeight','Demi' )
               ylabel('corr score-cell', 'FontSize',12, 'fontWeight','Demi' )
               set(gca,'FontSize',12 ,'fontWeight','Demi')
                
                                 %        xlim(a1, [ 1.2*min(F_1)   1.2*max(F_1)  ] )
                                 %        ylim(a1, 1.2*[ min(-PC_1) max(-PC_1)   ])
                                 %        ,'o','Color',[.7 .7 .7],'markerfacecolor', [.7 .7 .7] 
               
               
               
               
               
               
               
               
%              subplot(212)
%               plot( PC_a_sorted ,p_vals(ind),'.')
%               xlabel('loadings')
%               ylabel('pval corr score-cell  ')
%               ylim([0 1])
%               
     
        else 
            h_1 = 1;   
        end
        
      
        
end
