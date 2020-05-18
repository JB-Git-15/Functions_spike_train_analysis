function [DistV]=CorrelationVsShank(filename,tW,nJ,minRate,i_cov,shankOrder,MakePlot)

% This function calculates the Corr. or Cov. Matrix of the data in filename
% and splits its elements according to the distance between the shank in
% which in each neuron was recorded. Input variables tW,nJ,minRate (optional, 
% defaults to zero) are for computing the centered Spk Count Matrix using 
% SpkCountMat_Centered.m ; i_cov (optional, defaults to zero) specifies calculating
% the covariance, instead of the correlation matrix. shankOrder (optional, 
% defaults to the shanks listed in Ind) is a 1D array specifying shank order in case the
% shanks in Ind are not in the right order. MakePlot (optional, defaults 0)
% specifies if results should be plotted
% The function outputs DistV, a 1D cell array where each array member contains
% all the elements of the Corr/Cov matrix between pairs at a given
% distance, and CorrV, a 2D cell where each array element contains all the
% elements of the Corr/Cov matrix for a particular pair of shanks. 
% For instance, DistV{3} contains the CC of all the pairs recorded at a
% distance = 2 (e.g., shanks 2 and 4). CorrV{2}{4} contains the CC of all 
% the pairs in shanks 2 and 4.
%
% Example usage
% dV=CorrelationVsShank('SpkCells_Act6650_7450_NoStim',1,5,0,0,shankOrder,1);

load(filename)

% This loads a .mat file with our standard format: a cell array called Mat in 
% which each element is the spike train of 1 neuron (in s), and a 2D array called 
% Ind with the tetrode/shank and cluster of the corresponding cell.

% Set default input variables
if nargin<4
    minRate=0;
end

if nargin<5
    i_cov=0;
end

if nargin<6
    shankOrder=unique(Ind(:,1));
end

if nargin<7
    MakePlot=0;
end

% Calculate 'centered' spike count matrix
[~,SpkMC,~]=SpkCountMat_Centered(filename,tW,nJ,minRate);

nCells=length(SpkMC(1,:));
nShanks=length(unique(Ind(:,1)));
Shanks=unique(Ind(:,1));

% Check shankOrder is compatible with the number of shanks in the data
if length(shankOrder)~=nShanks
    disp('!!')
    disp('Number of shanks in shankOrder different from number of shanks in Ind array!!')
    disp('!!')
    return
end

% Find neurons recorded in each shank and rearrange Spike Count Matrix
SpkMCsort=zeros(size(SpkMC));
counter=0;
curri=counter+1;
for is=1:nShanks
    shank=shankOrder(is);
    icells{is}=find(Ind(:,1)==shank);
    counter=counter+length(find(Ind(:,1)==shank));
    DistV{is}=[];
    SpkMCsort(:,curri:counter)=SpkMC(:,icells{is});
    curri=counter+1;
end

% Calculate correlation/covariance matrices
if i_cov
    CMsort=cov(SpkMCsort);
    CM=cov(SpkMC);
else
    CMsort=corrcoef(SpkMCsort);
    CM=corrcoef(SpkMC);
end

% Array with all the (non repeated) elements of the corr. matrix.
CMV=CM(triu(ones(nCells,nCells),1)>0);
xh=min(CMV):(max(CMV)-min(CMV))/50:max(CMV);
[n,x]=hist(CMV,xh);
aux=cumsum(n)/sum(n);
alpha=0.005;
iL=find(aux<(1-alpha/2),1,'last');
iF=find(aux>(alpha/2),1,'first');
clim=[x(iF) x(iL)];
hlim=[min(CMV) max(CMV)];


% PCA
[coeff,latent]=pcacov(CM);
pc1=coeff(:,1); %pc1 loadings

% Separate elements of Corr Matrix by distance between shanks
for rows=1:nShanks
    for cols=rows:nShanks
        irow=icells{rows};
        icol=icells{cols};
        currMat=CM(irow,icol);
        CMpieces{rows}{cols}=currMat; %corr matrix for current shank pair
        if cols==rows
           lth=length(currMat(:,1));
           upm=triu(ones(lth,lth),1);          
           DistV{1}=[DistV{1};currMat(upm>0)];
           Hpieces{rows}{cols}=currMat(upm>0);
        else
           dist=cols-rows+1;
           DistV{dist}=[DistV{dist};currMat(1:numel(currMat))'];
           Hpieces{rows}{cols}=currMat(1:numel(currMat));
        end
    end
end

%Calculate mean and std of CCs as a function of shank distance
for i=1:nShanks
    meanCC(i)=mean(DistV{i});
    stdCC(i)=std(DistV{i});
end


% PLOT RESULTS
if MakePlot
    
    lSep=0.005;
    lPlot=(0.9-lSep*(nShanks-1))/nShanks;
    
    figure('position',[204   560   790   740]) %Corr Matrices
    for i=1:nShanks       
        for j=i:nShanks
            nPlot=(i-1)*nShanks+j;
            a(nPlot)=axes('position',[0.05+((j-1)*(lPlot+lSep)) 0.95-((i*lPlot)+((i-1)*lSep)) lPlot lPlot]);

            imagesc(CMpieces{i}{j});
            caxis(clim);
            set(gca,'xtick',[],'ytick',[])
            if i==j
                xlabel(['Shank ' num2str(i)])
            end
        end        
    end
    a(nPlot+1)=axes('position',[0.1 0.1 0.3 0.3]); %eigenvalues
    plot(latent,'bo','markerfacecolor','b','markersize',3),hold on
    plot([0 nCells],[1 1],'k:')
    xlim([0 nCells+1])
    ylabel('Fraction of Variance')
    xlabel('PC')
    a(nPlot+1)=axes('position',[0.24 0.24 0.15 0.15]); %pc1 loadings
    plot(sort(pc1),'ro','markerfacecolor','r','markersize',3),hold on
    plot([0 nCells],[0 0],'k:')
    xlim([1 nCells])
    ylabel('PC1 loading')
    xlabel('Neuron')
    
    
    
    figure('position',[404   560   790   740]) %Histograms
    for i=1:nShanks        
        for j=i:nShanks
            nPlot=(i-1)*nShanks+j;
            a(nPlot)=axes('position',[0.05+((j-1)*(lPlot+lSep)) 0.95-((i*lPlot)+((i-1)*lSep)) lPlot lPlot]);

            dh=(hlim(2)-hlim(1))/20;
            xh=hlim(1)+dh/2:dh:hlim(2)-dh/2;
            [n,nx]=hist(Hpieces{i}{j},xh);
            [x,y]=MakeHistLines(nx,n);
            plot(x,y),hold on
            yl=get(gca,'ylim');
            plot([0 0],yl,'r')
            xlim(clim);
            set(gca,'xtick',[],'ytick',[])
            if(i==j)
               set(gca,'xtick',[clim(1) 0 clim(2)]);
               xlabel(['CC (Shank ' num2str(i) ' )']);
            end        
        end      
    end
    
    figure('position',[604   560   790   740]) %Distance dependence
    for i=1:nShanks  
        nPlot=i;
        a(nPlot)=axes('position',[0.05+((i-1)*(lPlot+lSep)) 0.95-lPlot lPlot lPlot]); %Histograms for each distance
        
        [n,nx]=hist(DistV{i},xh);
        [x,y]=MakeHistLines(nx,n);
        plot(x,y),hold on
        yl=get(gca,'ylim');
        plot([0 0],yl,'r')
        xlim(clim);
        xlabel(['CC (Distance ' num2str(i-1) ' )'])
        set(gca,'ytick',[],'xtick',[0]);
        if (i==1 || i==nShanks)
            set(gca,'xtick',[clim(1) 0 clim(2)]);
        end        
    end 
   
    a(nPlot+1)=axes('position',[0.05 0.15+3*lPlot 0.9 2.5*lPlot]); % Mean and Std of histogram for each distance
    plot(0:nShanks-1,meanCC,'ko-','markerfacecolor','b','markersize',8),hold on
    plot(0:nShanks-1,stdCC,'ko-','markerfacecolor','r','markersize',8)
    plot([0 nShanks-1],[0 0],'k:')
    xlim([-0.5 7.5])
    minY=min([min(meanCC) min(stdCC)]);
    maxY=max([max(meanCC) max(stdCC)]);
    if minY<0
        ylim(1.25*[minY maxY]);
    else
        ylim([-0.25 1.25]*maxY);
    end
    set(gca,'xtick',[0:7])
    xlabel('Inter-Shank distance')
    ylabel('CC')
    legend('Mean CC','STD CC')
    
    a(nPlot+1)=axes('position',[0.05 0.1 0.9 3*lPlot]); % pc1 loadings for each shank
    for i=1:nShanks
        loadings=pc1(icells{i});
        Xshank=i+0.15*(rand(1,length(loadings))-0.5);
        for j=1:length(loadings) %silly trick to make dots opaque
            plot(Xshank(j),loadings(j),'ko','markerfacecolor','c','markersize',8),hold on
        end
        plot([0 nShanks+1],[0 0],'k:')
        xlim([0.5 nShanks+0.5])
        set(gca,'xtick',[1:nShanks])
        xlabel('Shank')
        ylabel('PC1 loading')
    end

end




