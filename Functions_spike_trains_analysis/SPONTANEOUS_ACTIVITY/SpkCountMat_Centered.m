function [SpkCountMat,SpkCountMatC,iLow] = SpkCountMat_Centered(filename,tWin,N_jitt,minRate)

% This function computes a 'centered' spike count matrix from spike
% train data using count windows of length tWin. Spike counts are centered 
% by the mean firing rate in windows of length tWin*N_jitt. A
% minimum Firing Rate threshold (minRate, optional -- defaults to 0) is allowed. Ouputs the centered
% (spkC) and conventional (spk) spike count matices and an array (iLow) of indices of the
% neurons excluded for not meeting the firing rate criterium.
%
% Example usage:
% [spk,spkC,iLow] = SpkCountMat_Centered('SpkCells_c037ActAll0_575',1,5,1);

load(filename)

% This loads .mat file with our standard format: a cell array called Mat in 
% which each element is the spike train of 1 neuron (in s), and a 2D array called 
% Ind with the tetrode/shank and cluster of the corresponding cell.

% Set default input variables
if nargin<4
    minRate=0;
end

nCell=length(Mat);

% Calculate time range where spikes are contained
maxt=zeros(1,nCell);
mint=zeros(1,nCell);
for i=1:nCell
    maxt(i)=max(Mat{i});
    mint(i)=min(Mat{i});
end
t_range(1)=0.001*round(1000*min(mint))-tWin;
t_range(2)=0.001*round(1000*max(maxt))+tWin;

% Calculate spike count matrix (1 column per neuron).
t_bin=t_range(1)+tWin/2:tWin:t_range(2)-tWin/2;
SpkCountMat=zeros(length(t_bin),nCell);

for i=1:nCell
    
    spk= Mat{i};
    spk= spk(  (spk >= t_range(1)  ) );
    spk= spk(  (spk <= t_range(2)  ) );   
    SpkCountMat(:,i) = hist(spk,t_bin);
end

% Make sure Spike count matrix is an integer number of multiples of N_jitter
nBins=size(SpkCountMat,1);
nBins_o=nBins;
nJits=floor(nBins/N_jitt);
if nBins_o>nJits*N_jitt
    SpkCountMat((nJits*N_jitt)+1:end,:)=[];
    nBins=size(SpkCountMat,1);
end

% Select only neurons with mean firing rates >= minRate
iLow=find((mean(SpkCountMat)/tWin)<minRate);
disp(' ')
if ~isempty(iLow)
    disp(['Cells ' num2str(iLow) ' were excluded because they fire'])
else
    disp(['No Cells were excluded because they fire'])    
end
disp(['less than ' num2str(minRate) ' Hz across the whole recording'])
disp(' ')
SpkCountMat(:,iLow)=[];
nCell = nCell-length(iLow);

% Local centering on Spike count matrix. Roughly equiv to high-pass.
SpkCountMatC=zeros(size(SpkCountMat));
nn=N_jitt;
    while nn<=nBins

        chunk=SpkCountMat(nn-N_jitt+1:nn,:);
        mChunk=mean(chunk);
        SpkCountMatC(nn-N_jitt+1:nn,:)=chunk-repmat(mChunk,N_jitt,1);

        nn=nn+N_jitt;

    end
 
end