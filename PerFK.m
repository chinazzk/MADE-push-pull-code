function [ PerFK ] = PerFK( K )
%This code transforms a matlab 3d k matrix into a column matrix for input to
%powerflow. Matlab matrix structure is LeftRight-TopBottom-pages(z), and powerflow matrix structure
%is LeftRight-frontback-layers(xy).Written by Hamed D. Ibrahim on Feb. 9, 2012.
%McCabe complexity is 2 [mlint('PerFK','cyc')]

K=exp(K); %fractal conditioning output is natural log (InK)
m=size(K,1);%Length of x axis
n=size(K,2);%Length of y axis
p=size(K,3);%Length of z axis
PerFK=zeros(m*n*p,1);%Pre allocation


for a=1:p%Length of z axis
    M=flipud(K(1:m,1:n,a));%Flip matrix upside down
    M=M';%Transpose to prepare matrix for vector(column)arrangement
    M=reshape(M,m*n,1);%Turn matrix into a column vector
     
        
  PerFK((((a-1)*m*n)+1):(((a-1)*m*n)+(m*n)))=M;%Replace zeros in PerFK 
                                               %with column vector
                                               %representing each z layer
  
        
end
        



end

