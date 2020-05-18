function [Lambda,Coefs,CosAngles]=PCA_TimeBin(filename,nJ,minRate,tW,MakePlot,iCells)

% This function explores the dependency of the PCA decomposition of a
% population recording during spontaneous activity on the timeBin used for
% computing spike counts. Inputs are the name of the input data file, and
% the number of bins nJ used for 'centering' the spike count vectors.
% Optional inputs are a minimum firing rate threshold (minRate, defaults to zero) for
% including in the population data set, the set tW of timeBins used for the
% analysis (1D array, defaults to 0.005 0.025 0.05 0.1 0.25 0.5 1 s), and a binary
% integer MakePlot that specifies if results should be plotted. The
% function outputs cell arrays Lambda and Coefs containing eigenvalue decomposition and the loadings of the
% first 9 PCs, for each timeBin. Also outputs CosAngles, a 3D array
% containing the cosine of the angles between each of the first 9 PCs (fist
% dimension) between pairs of timeBins (2nd and 3rd dimensions).
% 
% Example usage
% [L,C,A]=PCA_TimeBin('SpkCells_c037ActAll0_575',5,1,[0.01 0.1 1],1);

load(filename)

% This loads .mat file with our standard format: a cell array called Mat in 
% which each element is the spike train of 1 neuron (in s), and a 2D array called 
% Ind with the tetrode/shank and cluster of the corresponding cell.

% Set default input variables
if nargin<3
    minRate=0;
end

if nargin<4
    tW=[0.005 0.025 0.05 0.1 0.25 0.5 1];
end

if nargin<5
    MakePlot=0;
end

if nargin<6
    iCells=[];
end

% Fixed Parameters
nPC=9;

nT=length(tW);
Norm=zeros(nT,nPC);
for i=1:nT
    
    disp(['Processing timBin = ' num2str(tW(i)) ' s'])
    
    [~,SpkMC,~]=SpkCountMat_Centered(filename,tW(i),nJ,minRate,[],iCells);
    
    nCells=length(SpkMC(1,:));
    
    %standardize data
    stdData=std(SpkMC);
    SpkMCs=SpkMC./repmat(stdData,length(SpkMC(:,1)),1);
    
    [c,~,l]=pca(SpkMCs);
    
    Lambda{i}=l;
    Coefs{i}=c(:,1:nPC);
    Norm(i,:)=sum(c(:,1:nPC))/sqrt(nCells);
end

CosAngles=zeros(nPC,nT,nT);
for i=1:nT
    c1=Coefs{i};
    for j=i:nT
        c2=Coefs{j};
        for k=1:nPC
             CosAngles(k,i,j)=c1(:,k)'*c2(:,k); 
             CosAngles(k,j,i)=CosAngles(k,i,j); 
        end        
    end
end

if MakePlot % Plot results
    
    Gmin=0;
    Gmax=0.75;
    figure('position',[175   717   721   615])
    for i=1:nT
        l=Lambda{i};
        G(i)=Gmin+(i-1)*((Gmax-Gmin)/(nT-1));
        subplot(1,2,1)
        plot(l,'marker','o','color',G(i)*[1 1 1],'markerfacecolor',G(i)*[1 1 1]),hold on
        xlim([0 nCells])

        subplot(1,2,2)
        plot(l,'marker','o','color',G(i)*[1 1 1],'markerfacecolor',G(i)*[1 1 1]),xlim([0 10]),hold on
    end
   
    figure('position',[905   715   668   618],'name','Abs(Cos([angle])) for PCs for different tW')
    for i=1:nPC
        subplot(3,3,i)
        imagesc(squeeze(abs(CosAngles(i,:,:))));
        colormap('gray')
        caxis([0 1]);
        title(['PC ' num2str(i)])
        set(gca,'xtick',[],'ytick',[]);
    end
    
    figure('position',[1585 519 355 812],'name','Abs(Cos([angle])) with the uniform vector')
    for i=1:nPC
        yp=(0.9-((nPC-1)*0.01))/nPC;
        yc=0.95-i*yp-(i-1)*0.01;
        axes('position',[0.1 yc 0.8 yp])
        plot(tW,Norm(:,i),'ko-','markerfacecolor','k'),hold on
        plot(tW,zeros(1,length(tW)),'k:')
        ylim([-1 1])
        xlim([tW(1) tW(end)])
        set(gca,'ytick',[0 1])
        set(gca,'xtick',[0.01 0.1 1])
        set(gca,'xticklabel',[])
        ylabel(['PC ' num2str(i)])
        if i==nPC
            set(gca,'xticklabel',[0.01 0.1 1])
            xlabel('Spike Count Bin (s)')
        end
        set(gca,'xscale','log')
        
    end
    
    
    
end