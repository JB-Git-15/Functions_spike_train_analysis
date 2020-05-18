function  Plot_raw_data(filename )

%%%  This function allows to see the raw data in spikes
%%%  reads the Mat and trace alll the simultaneous activity of the network 
%%%  in windows of four seconds. 





        load(filename);
         spikes  = Mat;
%%%% Default parameters        
        Window_length_in_s = 10;
%        remain             =    assignopts (who, varargin);
        
        N_neurons          =    length(spikes);
 
% cherche min and max total...
        min_  = 10000;
        max_  = 0;
     for n_neur =  1 : N_neurons
         min_ = min([min_; spikes{1,n_neur}]);  
         max_ = max([max_; spikes{1,n_neur}]);
     end
        
     
    for time = min_ : Window_length_in_s  :  max_
        
        init_time  = time;
        end_time   = time + Window_length_in_s;
        
        scrsz = get(0, 'ScreenSize');
        
        n_neuron_i       =   ones(N_neurons,1);
        figure('Color','white', 'Position', [1 scrsz(4) scrsz(3) scrsz(4)] );
          
        
        
          hold on
              for n_neur =  1 : N_neurons
                  
 
                    spikes_in_wind = spikes{1,n_neur};
                    spikes_in_wind = spikes_in_wind(spikes_in_wind > init_time & spikes_in_wind <= end_time);
                    
                    if length(spikes_in_wind) 
                        for i = 1 : length(spikes_in_wind) 
                           plot([spikes_in_wind(i)  , spikes_in_wind(i)], [(n_neur - .9) n_neur  ] ,'-k' )
                        end
                    end
              end 
          hold off  
           xlabel('time (s)')
           ylabel('neurons')
           
           pause
           
           
           close
    
    end
        
     
  

