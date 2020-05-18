function [ left   bottom  width  height] = Dimensionate_frame(All_data,i,j)


s        = size(All_data);
n        = s(1);     % 
m        = s(2);
 
f_dy     = .29;
height   =  1/ (n*(1 + f_dy ) + 3*f_dy );   
dy       = height*f_dy;   

f_dx     = .15;

width    = 1/ (m*(1 + f_dx ) + 4*f_dx  ); 
dx       = width*f_dx;   

  

left   =    4*dx    + (j-1)*( width +dx );
bottom =    3*dy      + (n-i)*( height+dy );  

 

end

%   Key commands associated with this function :

%   dimS      = get(0,'ScreenSize');  
%   figure('Color','white','position',dimS)
%   a(nPlot+1)=axes('position',[ left   bottom  width  height]);
            
 









%   Original  : 
% s        = size(All_data);
% n        = s(1);     % 
% m        = s(2);
%  
% f_dy     = .29;
% height   =  1/ (n*(1 + f_dy ) + f_dy );   
% dy       = height*f_dy;   
% 
% f_dx     = .15;
% width    = 1/ (m*(1 + f_dx ) + f_dx/2  ); 
% dx       = width*f_dx;   
% 
%   
% 
% left   =    dx/2  + (j-1)*( width +dx );
% bottom =    dy    + (n-i)*( height+dy );  


