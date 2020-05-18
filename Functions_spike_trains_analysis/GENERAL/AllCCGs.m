function [ccgs,t]=AllCCGs(Mat,SR,included,epochs,bin,Nbin)

        %epochs = round(epochs*SR);
        nCells = length(Mat);
        iTot   = 1;
        rs     = [];
        cl     = [];
            for i=1:nCells
                        rs=[rs;round(Mat{i}*SR)];
                        cl=[cl;(i+1)*ones(length(Mat{i}),1)];
                if included(i)==1
                        gsub(iTot)=i+1;
                        iTot=iTot+1;
                end
            end
             
 [ccgs,t] = CCG(rs,cl,bin,Nbin,SR,gsub,'scale',epochs);    % 'scale'  ,'hz'
 
 
 
end