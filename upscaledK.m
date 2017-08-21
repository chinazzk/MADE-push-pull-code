function [ a ] = upscaledK(K,r,s,t)
%UNTITLED This code implements a rectangular cuboid scaling by computing sub
%block averages of a domain grid. 
%Cuboid dimensions are positive factors of domain grid lengths (the X,Y,Z
%lengths of your domain).
%a=original domain grid with block averaging.
%b=scaled domain grid.
%McCabe complexity is 4 [mlint('upscaledK','cyc')]



m=size(K,1);%Length of x axis
n=size(K,2);%Length of y axis
p=size(K,3);%Length of z axis

%A=[r s t];%This is the block size in the x, y, and z axis respectively.
            

Ktilde=zeros(m,n,p);%Pre-allocation of domain grid 
%Kscaled=zeros(m/r,n/s,p/t);%Pre-allocation of scaled domain

for j= 1:p/t;%z dimension of block
    
    for w=1:n/s;%y dimension of block
       
    for z= 1:m/r;%x dimension of block
           
        Kz= K((((z-1)*r)+1):(z*r), (((w-1)*s)+1):(w*s), (((j-1)*t)+1:(j*t)));
        Ktilde((((z-1)*r)+1):(z*r), (((w-1)*s)+1):(w*s), (((j-1)*t)+1:(j*t)))=mean2(Kz);
        %Kscaled(z,w,j)=mean2(Kz); 
        
    end
    
    end
    
      
end

a=Ktilde;
%b=Kscaled;%Replace grid with upscaled values
%save('Kscaled.txt','-ascii','Kscaled');%This saves a copy of the scaled hydraulic conductivity

  

end

