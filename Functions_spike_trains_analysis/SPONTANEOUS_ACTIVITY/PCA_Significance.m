function [SpkMCs,eigsData,eigsPermMat]=PCA_Significance(filename,tW,nJ,minRate,Nperm,iProjOut,MakePlot)

% This function calculates how many PCs are significant from a populatio recording of spontaneous activity.
% In addition to the input data (see below), the function needs tW, the
% time bin for calculating spike counts, and nJ, so that spike counts are
% centered on windows of length nJ*tW. Surrogate data will be generated by
% locally permuting nJ spike counts. minRate (optional, defaults to zero)
% is a minimum firing rate threshold for including cells into the spike
% count matrix. Nperm is the number of surrogate data sets, and iProjOut is
% a 1D array specifying the dimensions (any set of integers between 1 and
% nCells) to be 'projected out' before shuffling. MakePlot (defaults to zero) is a binary
% integer specifying whether results should be plotted. Our algorithm is a
% modification for spike trains of the procedure described in Schwartz et
% al., J. Vis. 2006.
%
% Example usage:
% [spk,e,eP]=PCA_Significance('SpkCells_c037ActAll0_575',1,5,1,100,[1:6],1);
% Produces the centered Spike Count Matrix, spk, the igenvalues, e, of the PCA
% decomposition of the data and a matrix eP with one eigenvalue
% decomposition for each of the Nperm=100 surrogates, where only the data
% projected onto the subspace orthogonal to the first 6 PCs is shuffled.
% Spike counts are made in windows of tW=1 sec, centered on 5 sec windows
% and only neurons with an overall firing of > 1Hz are considered.

load(filename)

% This loads .mat file with our standard format: a cell array called Mat in 
% which each element is the spike train of 1 neuron (in s), and a 2D array called 
% Ind with the tetrode/shank and cluster of the corresponding cell.

% Set default input variables
if nargin<4
    minRate=0;
end

if nargin<5
    Nperm=1;
end

if nargin<6
    iProjOut=[];
end

if nargin<7
    MakePlot=0;
end

% Calculate 'centered' spike count matrix
[~,SpkMC,~]  = SpkCountMat_Centered(filename,tW,nJ,minRate);

nCells       = length(SpkMC(1,:)); % number of neurons passing the minRate threshold
nJits        = length(SpkMC(:,1))/nJ; % number of 'isolated' time chunks

% Standardize data

stdData      = std(SpkMC);
SpkMCs       = SpkMC./repmat(stdData,length(SpkMC(:,1)),1);

%Do pca of standardized data
[cfData,scData,eigsData] = princomp(SpkMCs);

%Project data onto the space orthogonal to the pc's specified in iProjOut
iRetain           = 1:nCells;
iRetain(iProjOut) = [];
SpkMCsProj        = scData(:,iRetain)*cfData(:,iRetain)'; %this data matrix has dimensions nTimeBins*nCells, but has rank nCells - length(iProjOut)

% Generate and do PCA on surrogate data 
inds         = zeros(size(SpkMCs));
SpkMCsPerm   = zeros(size(SpkMCs));
eigsPerm     = zeros(size(eigsData));
eigsPermMat  = zeros(nCells,Nperm);
for i=1:Nperm
    % first generate all local random permutations
    allPerms = arrayfun(@(x)randperm(nJ),(1:nJits*nCells)','UniformOutput',0);
    allPermMat = cell2mat(allPerms)'; % this matrix is nJ rows x nJits*nCells columns
    % make matrix of permuted indices
    for j=1:nJits
        inds(1+(j-1)*nJ:j*nJ,:)=allPermMat(:,1+(j-1)*nCells:j*nCells)+((j-1)*nJ*ones(nJ,nCells));
    end
    % make permutated data matrix SpkMCsPerm: still locally centered and standardized
    for j=1:nCells
        SpkMCsPerm(:,j)=SpkMCsProj(inds(:,j),j);
    end  
    SpkMCsPProj=SpkMCsPerm*cfData; %Permutated data in the orig. data PC axis
    SpkMCsPProj(:,iProjOut)=[]; %Projection of the data on to the iRemain subspace
    [~,~,eigsPermPre]=princomp(SpkMCsPProj); %order variances within this subspace
    eigsPerm(iRetain)=eigsPermPre*(sum(eigsData(iRetain))/sum(eigsPermPre)); %scale so total variance in iRetain is the same as before shuffling
    eigsPerm(iProjOut)=NaN; %for plotting
    eigsPermMat(:,i)=eigsPerm;
end

if MakePlot
    figure('Color','white')
    hold on
    plot(eigsPermMat,'r')    
    
    plot(eigsData,'b.'),


    hold off
    xlabel('# PC', 'FontSize',18)
    ylabel(' % Explained variance', 'FontSize',18)
    set(gca, 'FontWeight','bold')
    
end
