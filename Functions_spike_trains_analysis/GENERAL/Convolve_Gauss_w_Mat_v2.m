function [rslt,time] = Convolve_Gauss_w_Mat_v2(Mat,sdt_1,sdtmax,dt,tlim)

% Function that convolves (smooths) a cell array Mat containing a point process on each
% elemment with a Gaussian kernel of width std sampled at resolution dt.
% The smoothed output is evaluated from tlim(1) until tlim(2) (at
% resolution dt). To avoid edge effects, the smoothed point processess are
% transiently extended beyond tlim (by an amount proportional to stdmax).


% In this new version the kernel is only computed once ! 

if nargin<5
    nCell = size(Mat,2);
    maxt=zeros(1,nCell);
    mint=zeros(1,nCell);
    for i=1:nCell
        maxt(i)=max(Mat{i});
        mint(i)=min(Mat{i});
    end
    tlim(1)=min(mint);
    tlim(2)=max(maxt);
end


nts    = round((tlim(2)-tlim(1))/dt);
tvec   = tlim(1)-5*sdtmax:((tlim(2)-tlim(1))/nts):tlim(2)+5*sdtmax;
nCells = length(Mat);
csign  = zeros(nCells,length(tvec));

N_time_steps           = length(tvec);
time_centered          = tvec - tvec(1);

time_centered_extended     = sort(-time_centered);
time_centered_extended     = time_centered_extended(1:end-1);
time_centered_extended     = [time_centered_extended,time_centered]; % 2*N_time_steps - 1 elements
extended_impulse_response  = normpdf(time_centered_extended,0,sdt_1);

  parfor i=1:nCells

    spkt=Mat{i}(Mat{i}>tvec(1) & Mat{i}<=tvec(end));
    s=size(spkt);
    
    if s(1)==1 && s(2)>1
        spkt=spkt';
    end
    
    if length(spkt)>1
 
        buffer = zeros(size(tvec)); 
        for n_spike_t = 1 : length(spkt)
             mu              =  spkt(n_spike_t);
             [val,ind]       =  min( (mu - tvec).*(mu - tvec));
              
              Impulse_reponse=  extended_impulse_response(N_time_steps- ind + 1 :2*N_time_steps- ind);             
             % Impulse_reponse =  normpdf(tvec,mu,sdt_1);     
             buffer          =  buffer + Impulse_reponse; 
        end
        csign(i,:) =  buffer;           
        
    elseif length(spkt)==1
        csign(i,:) = normpdf(tvec,spkt,sdt_1)   ;          
    else
        csign(i,:) = zeros(1,length(tvec));
    end
    
end

rslt = csign(:,tvec>tlim(1) & tvec<=tlim(2))';
time = tvec(tvec>tlim(1) & tvec<=tlim(2));


end
