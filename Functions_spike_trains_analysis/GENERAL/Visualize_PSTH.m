function Visualize_PSTH(PSTH,time,varargin)



      rank          =    'none';     % time before stimulus ( en s)
      remain        =    assignopts (who, varargin);

    
      
 if strcmp(rank,'fr')  % sort by firing rate
      
   [val,N_n] =  sort(sum(PSTH,2),'descend');
    PSTH     =  PSTH(N_n,:);
    
end


   figure('Color','white')
      imagesc(time, 1 : size(PSTH,1), PSTH)
      xlabel('time (s)')
      ylabel('neurons')
      title(' averaged mean cell response to stimulus (Hz)')
      colorbar





end
      