function [ S ] = KrigsFast3d(pts,cov,m)
% Kriging script by NMonnig 6/8/2012
% Interpolates using simple kriging using all conditioning points
% shifted to zero mean
% WARNING:  positive vs negative directionality has not been verified
%           OK if using a symmetric covariance structure
%           else, check for consistency in "upstream vs downstream"

% make it zero mean for simple kriging
mutemp=mean(pts(:,4)); pts(:,4)=pts(:,4)-mutemp; 
midx=m(1);midy=m(2);midz=m(3);
np=length(pts(:,1));

% C = matrix of covariances between conditioning points
C=zeros(np,np);
for i=1:np
    for j=1:np
        C(i,j)=cov(pts(i,1)-pts(j,1)+1+midx,pts(i,2)-pts(j,2)+1+midy,pts(i,3)-pts(j,3)+1+midz);
    end
end

% Solve linear system to find weights
w=C\pts(:,4);

S=zeros(size(cov));
% Interpolate new points given weights by adding shifted cov array
for l=1:np
    S=S+w(l)*circshift(cov,[pts(l,1)-1 pts(l,2)-1 pts(l,3)-1]);
end

S=S(midx+1:2*midx,midy+1:2*midy,midz+1:2*midz)+mutemp;

