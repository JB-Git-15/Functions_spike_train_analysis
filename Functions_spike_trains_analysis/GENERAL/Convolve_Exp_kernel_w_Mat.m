function [rslt,time] = Convolve_Exp_kernel_w_Mat(Mat,tau,tau_max,dt,tlim)

% Function that convolves (smooths) a cell array Mat containing a point process on each
% elemment with an exponential kernel of tau sampled at resolution dt.
% The smoothed output is evaluated from tlim(1) until tlim(2) (at
% resolution dt). To avoid edge effects, the smoothed point processess are
% transiently extended beyond tlim (by an amount proportional to tau_max).

nts          = round((tlim(2)-tlim(1))/dt);
tvec         = tlim(1)-5*tau_max:((tlim(2)-tlim(1))/nts):tlim(2)+5*tau_max;
nCells       = length(Mat);
N_time_steps = length(tvec); 
csign        = zeros(nCells,length(tvec));

%%% In this version we calculate only once the kernel...
time_temp            = tvec - tvec(1);
Impulse_reponse_max  = tau*exp(-time_temp/tau);


 for i=1:nCells

    spkt=Mat{i}(Mat{i}>tvec(1) & Mat{i}<=tvec(end));
    s=size(spkt);
    
    if s(1)==1 && s(2)>1
        spkt=spkt';
    end
    
    if length(spkt)>1
 
        buffer = zeros(1,N_time_steps ); 
        for n_spike_t = 1 : length(spkt)
             t_s               =  spkt(n_spike_t);           % time spike
             [val,Index_closest_spk] =  min( (tvec - t_s).*(tvec - t_s) ); % index of the closest spike element 
             
             Impulse_reponse                                    =  zeros(1,N_time_steps); 
             Impulse_reponse(Index_closest_spk : N_time_steps)  =  Impulse_reponse_max(1 : N_time_steps- Index_closest_spk+1)  ;     
             
             buffer          =  buffer + Impulse_reponse; 
        end
        csign(i,:) =  buffer;           
        
    elseif length(spkt)==1
             buffer = zeros(1,N_time_steps ); 
            for n_spike_t = 1 : length(spkt)
             t_s               =  spkt(n_spike_t);           % time spike
             [val,Index_closest_spk] =  min( (tvec - t_s).*(tvec - t_s) ); % index of the closest spike element 
             
             Impulse_reponse                                    =  zeros(1,N_time_steps); 
             Impulse_reponse(Index_closest_spk : N_time_steps)  =  Impulse_reponse_max(1 : N_time_steps- Index_closest_spk+1)  ;     
             buffer                                             =  buffer + Impulse_reponse;
            end
             
             csign(i,:)                                         =  buffer;
    else
        csign(i,:) = zeros(1,length(tvec));
    end
    
end

rslt  =  csign(:,tvec>tlim(1) & tvec<=tlim(2))';
time  =  tvec(tvec>tlim(1) & tvec<=tlim(2));


end
