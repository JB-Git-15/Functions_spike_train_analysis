function [t_1 t_2] = Determine_time_limits(fime_name)

%%%%% Attention ! changer le Mat par spikes ou l'inverse : (le standard va être Mat)
load(file_name)

N_neurons = length(Mat);

                 mini = 100000;
                 maxi = -10;
                 
            for i = 1 : N_neurons
                 mini  = min(mini, Mat{1,i}(1));
                 maxi  = max(maxi, Mat{1,i}(end));    
            end

            
            t_1 = mini;
            t_2 = maxi;
            
end