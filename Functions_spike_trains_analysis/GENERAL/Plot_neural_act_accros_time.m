function h_1 =   Plot_neural_act_accros_time(filename,Gaussian_kernel_std,t_range,rateTh,List_to_exclude )
 
load(filename)
 
    N_neurons          = length(Mat);
   % rateTh             = .5; % min FRate threshold

    sdtmax             = 2*Gaussian_kernel_std;
    dt                 = Gaussian_kernel_std/4;
 
%   Calculate mean firing rate for each cell 

     Mean_FR  = zeros( N_neurons,1);
         
    for i= 1 : N_neurons
                       spk                         = Mat{i};
                       Mean_FR(i)                  = length(spk)/(t_range(2) -t_range(1));
    end
    
    [rslt,time] = Convolve_Gauss_w_Mat_v2(Mat,Gaussian_kernel_std,sdtmax,dt,t_range);
     rslt       = rslt';
   
    iLow         = find( Mean_FR < rateTh);
    
    
    
    %  temp         = rslt(iLow,:);
    %  rslt(iLow,:) = [];
             
   
    
    Number_dezens = ceil(N_neurons/10) ;%    
    a             = Number_dezens;
    b             = 10;
  
    
     %   N_neurons    = N_neurons - length(iLow);
 
        disp([ filename, '   '] )
        disp( ['Cells #' num2str(iLow') ' were excluded because they fire'] )
        disp( ['less than ' num2str(rateTh) ' across the whole recording'] )
        disp(' ')
        
  
         
        
        
        
        
 Dim_screen    =  get( 0, 'ScreenSize');
 Dim_screen(4) =  Dim_screen(4) - 120;
  
h_1 =  figure('position',Dim_screen,'name',['rateof units accross time : ', filename],'Color','w')  ; 
       
          count     = 1; 
           
         for i = 1 : a
             for j = 1 : b
                   
                 if count < ( N_neurons +1) 
                 
                 [ left  bottom width height] = Dimensionate_frame(ones(a,b), i,j) ;
                 
                
                    
                     h_temp         = axes('position',[ left   bottom   width  height]);
                     y_max          = 32;
                     
                     
                     if length(List_to_exclude)
                                 if ~prod(List_to_exclude - count)
                                                      plot(h_temp, time, rslt(count,:),'Color','r');

                                 else                 plot(h_temp, time, rslt(count,:));              

                                 end        
                     else     plot(h_temp, time, rslt(count,:))
                     
                     end
                                    % xlabel(' time (s)')
                                    % ylabel(['      ', num2str(count) ] )
                                      axis([t_range(1) t_range(2) 0 y_max ]);
                                       axis off
                     
                                      
                 end   

                  count = count + 1 ; 
             end
         end
        
        
        
        
    
end