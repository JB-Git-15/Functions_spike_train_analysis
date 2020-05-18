function Coh=LFP_Shank_Coherence_v2(file,nChan,SR,t_range,win,nFFT,freqs)

% Calculates coherence between channels in the lfp matrix. The lfp matrix
% is calculated reading from the raw (.eeg or .dat) data file. If you want to use
% an already created lfp matrix, use the _v1 version of this function.
% file is the file name. 
% nChan is the number of channels
% SR is the sampling rate in Hz
% t_range is a 2 element array with the initial and final times to be
% analyzed (in sec)
% win and nFFT are parameters of the mscohere matlab coherence function
% (check matlab help).
% freqs is a 2 element array with the minimum and maximum frequencies in
% the coherence plot (in Hz)
% The funciton plots matrix of average (across frequencies) coherences in the frequency range
% specified by freqs=[fmin fmax]

iSec=t_range(1); %initial time to be analyzed (in sec)
fSec=t_range(2); %final time to be analyzed (in sec)

 
arg1=[nChan round(SR*(fSec-iSec))];

lfp=bload(file,arg1,round(iSec*SR)*2*SR);

s=size(lfp);
nTraces=s(1);

for i=1:nTraces-1
    Coh(i,i)=1;
    for j=i+1:nTraces
        [c,f]=mscohere(lfp(i,:),lfp(j,:),win,[],nFFT,SR);
        ifr=find(f>=freqs(1) & f<=freqs(2));
        Coh(i,j)=mean(c(ifr));
        Coh(j,i)=Coh(i,j);
    end
end
Coh(nTraces,nTraces)=1;   

figure
imagesc(Coh)