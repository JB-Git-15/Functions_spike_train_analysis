function [deleteCells,statCells]=Examine_Stationarity(filename,MakePlot)

%This function analyzes the stationarity of the firing rate of each neuron
%in the population during the whole duration of the recording. Stationarity
%is assessed by looking at the CV of the spike count across time calculated
%with bins of 10 seconds. It outputs two 4-element cell arrays with the identity of the
%neurons with CV less than 0.25, 0.5, 0.75 and 1 (statCells) and its
%complement (the neurons that don't meet the condition). Optional argument MakePlot
% (default value 0) specifies if results should be plotted. The function plots
%the number of cells who meet each of the 4 conditions, the histogram of
%the rates of the cells who meet each of the 4 conditions, and the firing
%rate across time of each neuron in the data-set with a color corresponding
%to its CV.

load(filename)

% This loads .mat file with our standard format: a cell array called Mat in 
% which each element is the spike train of 1 neuron (in s), and a 2D array called 
% Ind with the tetrode/shank and cluster of the corresponding cell.

nCell=length(Mat);

%optional arguments
if nargin<2
    MakePlot=0;
end

% fixed parameters
tWin=10;
CVs=[0.25 0.5 0.75 1];

% make time_bin array for Spike count matrix
maxt=zeros(1,nCell);
mint=zeros(1,nCell);
for i=1:nCell
    maxt(i)=max(Mat{i});
    mint(i)=min(Mat{i});
end
tRange(1)=min(mint);
tRange(2)=max(maxt);

numChunks=floor((tRange(2)-tRange(1))/tWin);
t_bin=(tRange(1)-(tWin/2))+tWin*(1:numChunks);
t_end=t_bin(end)+tWin/2;
nBins=length(t_bin);

% Calculate spike count matrix 
SpkCountMat=zeros(nBins,nCell);
for i=1:nCell
    spk=Mat{i};
    spk(spk>t_end)=[];
    SpkCountMat(:,i)=hist(spk,t_bin);
end

% CV of spike count on windows of tWin=10 s
meanRates=mean(SpkCountMat);
cv=std(SpkCountMat)./meanRates;

% find cells with CV less than a few thresholds
for i=1:length(CVs)
    statCells{i}=find(cv<CVs(i));
    nStatCells(i)=length(find(cv<CVs(i)));
    allCells=1:length(cv);
    allCells(cv<CVs(i))=[];
    deleteCells{i}=allCells;
end



if MakePlot
    
    %classify cells according to CV
    colors=zeros(nCell,3);
    iblack=find(cv<CVs(1));
    iblue= cv>=CVs(1) & cv<CVs(2);
    colors(iblue,:)=repmat([0 0 1],length(find(iblue==1)),1);
    ired= cv>=CVs(2) & cv<CVs(3);
    colors(ired,:)=repmat([1 0 0],length(find(ired==1)),1);
    igreen= cv>=CVs(3) & cv<CVs(4);
    colors(igreen,:)=repmat([0 1 0],length(find(igreen==1)),1);
    icyan= cv>=CVs(4);
    colors(icyan,:)=repmat([0 1 1],length(find(icyan==1)),1);
    
    
    figure('position',[636         990        1520         340])
    subplot(1,5,1)
    plot(CVs,nStatCells,'bo'); %number of cells with CV less than x-axis
    xlim([0 1.2])
    xlabel('CV')
    ylabel('N Cells')
    subplot(1,5,2) %histograms of firing rate for each CV
    hist(meanRates(statCells{1})/tWin)
    xlabel('Firing Rate (Hz)')
    ylabel('N Cells')
    title(['CV = ' num2str(CVs(1))])
    subplot(1,5,3)
    hist(meanRates(statCells{2})/tWin)
    xlabel('Firing Rate (Hz)')
    ylabel('N Cells')
    title(['CV = ' num2str(CVs(2))])
    subplot(1,5,4)
    hist(meanRates(statCells{3})/tWin)
    xlabel('Firing Rate (Hz)')
    ylabel('N Cells')
    title(['CV = ' num2str(CVs(3))])
    subplot(1,5,5)
    hist(meanRates(statCells{4})/tWin)
    xlabel('Firing Rate (Hz)')
    ylabel('N Cells')
    title(['CV = ' num2str(CVs(4))])

    %Plot firing rates of each cell
    nSq=ceil(sqrt(nCell));
    lSep=0.01;
    lPlot=(0.9-lSep*(nSq-1))/nSq;
    xl=([tRange(1) t_end]);

    str=['Black CV<0.25 ; Blue 0.25<CV<0.5 ; Red 0.5<CV<0.75 ; Green 0.75<CV<1 ; Cyan CV>1'];
    figure('position',[367         530        1900         703],'name',str)
    for i=1:nSq

        for j=1:nSq
            nPlot=(i-1)*nSq+j;
            if nPlot<=nCell
                a(nPlot)=axes('position',[0.05+((j-1)*(lPlot+lSep)) 0.95-((i*lPlot)+((i-1)*lSep)) lPlot lPlot]);

                [x,y]=MakeHistLines(t_bin,SpkCountMat(:,nPlot));
                plot(x,y,'color',colors(nPlot,:))
                xlim(xl)
                xlabel(' ');ylabel(' ');set(gca,'yticklabel',[]);set(gca,'xticklabel',[]);
            else
                return
            end
        end

    end
end
