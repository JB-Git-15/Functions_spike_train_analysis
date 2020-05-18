function [ccgs,tccg] = CCGsFromMat(filename,bin,Nbin,SR,gsub,epochs,Normalisation,MakePlot)

% Function using Ken's CCG mex function for calculating CCGs fast.
%
% filename = input data (see below).
% bin = bin size in input units 
% Nbin = |max Lag| of the CCG divided by bin
% SR = sampling rate (optional). How many input units per second. Defaults
% to 1 (i.e., input units are sec)
% gsub = array with indices of neurons to be used (optional). Defaults to
% all the neurons in the array Mat.
% epochs = subset of times to be used for the CCG calculation (optional). 
% Defaults to time range spanned by the spikes. 
% MakePlot (optional, defaults to false)= binary integer specifying whether 
% results should be plotted.
%
% Example: If the spike times are in sec, then for making CCG with a bin of
% 2 ms and MaxLag of 100 ms, then do
% CCGsFromMat('SpkCells_Act6650_7450_NoStim',0.002,50);

load(filename)

% This loads .mat file with our standard format: a cell array called Mat in 
% which each element is the spike train of 1 neuron (in s), and a 2D array called 
% Ind with the tetrode/shank and cluster of the corresponding cell.

nCells=length(Mat);

if nargin<4
    SR=1;
end

if nargin<5 || isempty(gsub)
    gsub=1:nCells;
end

if nargin<6
    epochs=[];
end

if nargin<7
    MakePlot=0;
end
    
epochs=round(epochs*SR);

rs=[];
cl=[];
for i=1:nCells
    rs=[rs;Mat{i}*SR];
    %rs=[rs;round(Mat{i}*SR)]; %doesn't look like Ken's function needs integer times
    cl=[cl;i*ones(length(Mat{i}),1)];
end

%By default the normalization is by the product of the firing rates
[ccgs,tccg]=CCG(rs,cl,bin,Nbin,SR,gsub, Normalisation ,epochs);     % or CCG_Ken
  
if MakePlot % Plot matrix of plots with individual CCGs
    nSq=length(gsub);
    lSep=0.01;
    lPlot=(0.9-lSep*(nSq-1))/nSq;
        
    figure('position',[78   104   840   702])
    for i=1:nSq
        
        for j=i:nSq
        nPlot=(i-1)*nSq+j;
        a(nPlot)=axes('position',[0.05+((j-1)*(lPlot+lSep)) 0.95-((i*lPlot)+((i-1)*lSep)) lPlot lPlot]);
        
        [x,y]=MakeHistLines(tccg,ccgs(:,i,j));
        xl=[min(tccg) max(tccg)];
        plot(xl,[1 1],'r'),hold on
        plot([0 0],[0 1],'r')
        plot(x,y)
        xlim(xl)
        xlabel(' ');ylabel(' ');set(gca,'ytick',[]);set(gca,'xtick',[]);
        if(i==j)
           set(gca,'xtick',[xl(1) 0 xl(2)]);
           xlabel('lag');
        end
        
        end
        
    end
    
end
