function [x,y]=MakeHistLines(xb,yb)
%
%This function takes the xb and yb values of a histogram, to be 
%plotted with bar(xb,yb), and returns x,y so that plot(x,y) also gives 
%a bar plot (but not filled)
%
delx=xb(2)-xb(1);
x(1)=xb(1)-delx/2;
y(1)=0;
for i=1:length(yb)
    y(2*i)=yb(i);
    y(1+2*i)=yb(i);
    x(2*i)=xb(i)-delx/2;
    x(1+2*i)=xb(i)+delx/2;
end
x(2*(i+1))=xb(length(xb))+delx/2;
y(2*(i+1))=0;