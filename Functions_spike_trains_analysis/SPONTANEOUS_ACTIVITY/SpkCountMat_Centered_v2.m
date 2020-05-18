function [SpkCountMat,SpkCountMatC,iDelete,t_bin]=SpkCountMat_Centered_v2(filename,tWin,N_jitt,minRate,tRange,iCells)

% This function computes a 'centered' spike count matrix from spike
% train data using count windows of length tWin. Spike counts are centered 
% by the mean firing rate in windows of length tWin*N_jitt. A
% minimum Firing Rate threshold (minRate, optional -- defaults to 0) is allowed. 
% tRange (optional, defaults to min and max spike times in data set) specifies
% the initial and final times (in sec) to be considered for analysis for all neurons.
% iCells (optional, defaults to all cells in Mat) specifies the neurons
% that sould be used for analysis. Ouputs the centered (spkC) and conventional (spk) 
% spike count matices, an array (iDelete) of indices of the neurons excluded 
% and t_bin, the vector of time bins used for calculating the spike count
% histograms.
%
% Example usage:
% [spk,spkC,iDelete]=SpkCountMat_Centered('SpkCells_c037ActAll0_575',1,5,0.5,[],1:25);
% Will produce spike count matrices for neurons 1 to 25 in windows of 1 sec
% centered on 5 sec windows.

load(filename)

% This loads .mat file with our standard format: a cell array called Mat in 
% which each element is the spike train of 1 neuron (in s), and a 2D array called 
% Ind with the tetrode/shank and cluster of the corresponding cell.

nCell=length(Mat);

% Set default input variables
if nargin<4
    minRate=0;
end

if nargin<5 || isempty(tRange) %time range defined by min and max spike times in dataset
    maxt=zeros(1,nCell);
    mint=zeros(1,nCell);
    for i=1:nCell
        maxt(i)=max(Mat{i});
        mint(i)=min(Mat{i});
    end
    tRange(1)=min(mint);
    tRange(2)=max(maxt);
end

if nargin<6
    iCells=[];
end

% define time bin array. Number of bins should be an integer multiple of jitter windows
numChunks=floor((tRange(2)-tRange(1))/(tWin*N_jitt));
t_bin=(tRange(1)-(tWin/2))+tWin*(1:numChunks*N_jitt);
t_end=t_bin(end)+tWin/2;
nBins=length(t_bin);

% Calculate spike count matrix (1 column per neuron).
SpkCountMat=zeros(nBins,nCell);
for i=1:nCell
    spk=Mat{i};
    spk(spk>t_end)=[];
    SpkCountMat(:,i)=hist(spk,t_bin);
end

% Select neurons
if isempty(iCells) % population not specified in function call. Check minRate criterion
    iDelete=find((mean(SpkCountMat)/tWin)<minRate);
    disp(' ')
    if ~isempty(iDelete)
        disp(['Cells ' num2str(iDelete) ' were excluded because they fire'])
    else
        disp(['No Cells were excluded because they fire'])    
    end
    disp(['less than ' num2str(minRate) ' Hz across the whole recording'])
    disp(' ')
else  % use only cells specified in function call
    iDelete=1:length(SpkCountMat(1,:));
    iDelete(iCells)=[];    
end
SpkCountMat(:,iDelete)=[];  %remove cells not matching criterion
nCell=nCell-length(iDelete); %update number of cells

% Local centering on Spike count matrix. Roughly equiv to band-pass.
SpkCountMatC=zeros(size(SpkCountMat));
nn=N_jitt;
while nn<=nBins

    chunk=SpkCountMat(nn-N_jitt+1:nn,:);
    mChunk=mean(chunk);
    SpkCountMatC(nn-N_jitt+1:nn,:)=chunk-repmat(mChunk,N_jitt,1);
    
    nn=nn+N_jitt;

end
 
    
