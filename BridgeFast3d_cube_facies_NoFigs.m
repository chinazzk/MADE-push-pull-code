% Modified by Nate Monnig to allow different statistical properties for
% different facies 6/18/2012.
%
% Modified by David A. Benson to allow conditioning and 3-dimensions
% 11/25/2011.  Uses a companion function trend3d.m similar to kriging.
% Corrected feller trend surface to make zero mean surfacee which is simple kriging.
% Indexing is somewhat dyslexic.  Matrices(i,j,k) refer to x1, x2, x3 dirs.
% while matlabs rows are first position but "feel" like "y".  This affects
% plotting.  Also now allows different # nodes and delta x in all dirs. 

% Edited by Mark M Meerschaert on 14 April 2011
% Based on code by Hans-Peter Scheffler circa 2008
%
% Simulation code for fractal K field with operator scaling
% Code can produce a normal or stable (heavy tailed) log K field
% Field has stationary increments, like fractional Brownian motion
% Hurst index in x1 direction is H/a_1
% Hurst index in x2 direction is H/a_2
% IMPORTANT NOTE: The code requires 0 < H < a_i !!! 
%
% Input variable definitions
% alpha = 2 for normal noise, 0<alpha<2 for stable noise (heavy tail)
% H = Hurst index applied to filter
% A = Scale for noise variables (standard deviation if alpha=2)
% m = size of OUTPUT field (internal fields are 2M+1 by 2M+1)
% a1 = power law index of filter for x1 coordinate
% a2 = power law index of filter for x2 coordinate
% v1 = direction of x1 coordinate axes in radians 
% v2 = direction of x2 coordinate axes in radians 
% C1 = correlation length for x1 coordinate 
% C2 = correlation length for x1 coordinate  
% rho = 2 (see formula for filter phi)
% mu = sample standard deviation of log K field
% sigma = sample standard deviation of log K field
% nuse  = # of conditioning points to use for each generated point
% nreals = # of realizations to create
clear all
%matlabpool close
%matlabpool open % go parallel with default params

 % angles (azimuth, dip) associated with the scaling eigenvectors
 angle1=[0*pi/16,0*pi/16]; angle2=[8*pi/16,0]; angle3=[0*pi/16,8*pi/16];

 % Geometry of the output field (with current coordinates, center at (0,0))
 % maximum depth = 12.19 (depth to aquitard)
 origin=[-25,-25,3.6]; maximum=[25,25,12.19]; m = [128,128,512]; dx=(maximum-origin)./m; nreals=10;
 disp(strcat({'discretization blocks have size = '},num2str(dx(1)),{' by '},num2str(dx(2)),{' by '},num2str(dx(3))))

 %unit scaling eigenvectors in the i,j,k directions
 theta1=[cos(angle1(2))*cos(angle1(1)),cos(angle1(2))*sin(angle1(1)),sin(angle1(2))]; 
 theta2=[cos(angle2(2))*cos(angle2(1)),cos(angle2(2))*sin(angle2(1)),sin(angle2(2))]; 
 theta3=[cos(angle3(2))*cos(angle3(1)),cos(angle3(2))*sin(angle3(1)),sin(angle3(2))]; 
 
 % check to see if there is conditioning data
 if exist('conditioning3d.txt')==2
  docondition=1;
  cpts1 = load('conditioning3d.txt','-ascii');
  cpts1(:,1)=ceil((cpts1(:,1)-origin(1))/dx(1));  
  cpts1(:,2)=ceil((cpts1(:,2)-origin(2))/dx(2));  cpts1(:,3)=ceil((cpts1(:,3)-origin(3))/dx(3));
  cpts1=sortrows(cpts1, [1 2 3]);  cpts=cpts1;
  % consolidate repeated points in a single cell by averaging the conditioning points
  ncp=length(cpts(:,1)); k=ncp;
  while k>1
    loc=cpts(k,1:3);
    count=1;
    value=cpts(k,4);
    while k-count>0 && cpts(k-count,1)==loc(1) && cpts(k-count,2)==loc(2) && cpts(k-count,3)==loc(3)
        value=value+cpts(k-count,4);
        count=count+1;
    end
    cpts(k,4)=value/count;
    if count>1
        cpts(k-count+1:k-1,:)=[];
    end
    k=k-count;
  end
  ncp=length(cpts(:,1));  % final number of conditioning points
  disp(strcat({'number of conditioning points = '},num2str(ncp)));
 else docondition=0;  % no conditioning data
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp('Here we go...')
  tic;

% load facies data
load('3Dset.mat');
%column 1: y coordinates - towards north (if you want, you may add 166
%meter, to get to actual GPR coordinate system.)
%column 2: x coordinates - towards east (if you want, you may add 95
%meter, to get to actual GPR coordinate system.)
%column 3: depths for water table
%column 4: depths for water lower boundary of facies 1
%column 5: depths for lower boundary of facies 2
%column 6: depths for lower boundary of facies 3
%column 7: depths for aquitard (Eutaw clay)
ycoords=data(:,1);
xcoords=data(:,2);
wtable=data(:,3);
fac1bot=data(:,4);
fac3bot=data(:,5); %%%%%% ORDER SWITCHED IN DATA FILE
fac2bot=data(:,6);
fac4bot=data(:,7);
%fac5bot=data(:,8);


if abs(min(wtable)-max(wtable))>10^(-12)
    disp('water table is not flat!!!')
end
facies.origin=[min(xcoords),min(ycoords)]; 
facies.maximum=[max(xcoords),max(ycoords)]; 
facies.dx=[0.5,0.5]; % check this manually
facies.m =round((facies.maximum-facies.origin)./facies.dx)+[1 1]; 
facies.n=5; % model extends to top of aquitard (bottom of facies #5)

fac1bot=reshape(fac1bot,facies.m(1),facies.m(2));
fac2bot=reshape(fac2bot,facies.m(1),facies.m(2));
fac3bot=reshape(fac3bot,facies.m(1),facies.m(2));
fac4bot=reshape(fac4bot,facies.m(1),facies.m(2));
%fac5bot=reshape(fac5bot,facies.m(1),facies.m(2));

xcoords2=reshape(xcoords,size(fac1bot));
ycoords2=reshape(ycoords,size(fac1bot));
figure(21)
scatter3(origin(1)+dx(1)*(cpts(:,1)-1/2),origin(2)+dx(2)*(cpts(:,2)-1/2),...
    origin(3)+dx(3)*(cpts(:,3)-1/2),[],cpts(:,4))
hold on
colorshade=mean(cpts(:,4))*ones(size(xcoords2));
surf(xcoords2,ycoords2,fac1bot,colorshade)
%surf(xcoords2,ycoords2,fac2bot,colorshade)
%surf(xcoords2,ycoords2,fac3bot,colorshade)
%surf(xcoords2,ycoords2,fac4bot,colorshade)
%%surf(xcoords2,ycoords2,fac5bot,colorshade)
%hold off

 % open mat file with all facies data
 %facies=open('FACIES.mat');
 K_facies=cell(facies.n,1);
 %origin_facies=zeros(facies.n,3);
 
%%%%%%% determine which facies each point in the grid is in
 ID=zeros(m(1),m(2),m(3));
 for i=1:m(1)
    for j=1:m(2)
        % consider location of each block in ID array to be center of block
        ind=round( (origin(1:2)+[i-1/2 j-1/2].*dx(1:2)-facies.origin)./facies.dx)+[1 1];
        if ind(1)>0 && ind(1)<facies.m(1)+1 && ind(2)>0 && ind(2)<facies.m(2)+1  
            ID(i,j,:)=5;
            for k=m(3):-1:1
                z=origin(3)+dx(3)*(k-1/2);
                
                %if z <= fac5bot(ind(1),ind(2))
                    %ID(i,j,1:k)=5;
                    if z <= fac4bot(ind(1),ind(2))
                        ID(i,j,1:k)=4;
                        if z <= fac3bot(ind(1),ind(2))
                            ID(i,j,1:k)=3;
                            if z <= fac2bot(ind(1),ind(2))
                                ID(i,j,1:k)=2;
                                if z <= fac1bot(ind(1),ind(2))
                                    ID(i,j,1:k)=1;
                                end
                            end
                        end
                    end
                %end
                
            end
        end
    end
 end

% extend facies data to the rest of the cube 
first=find(ID,1,'first');
last=find(ID,1,'last');
[If Jf Kf]=ind2sub([m(1) m(2) m(3)],first);
[Il Jl Kl]=ind2sub([m(1) m(2) m(3)],last);
nonzero=ID(If:Il,Jf:Jl,Kf:Kl);
prepad=[If Jf Kf]-[1 1 1];
postpad=[m(1) m(2) m(3)]-[Il Jl Kl];
nonzero=padarray(nonzero,prepad,'replicate','pre');
ID=padarray(nonzero,postpad,'replicate','post');
%%%%%%%%%%%%%%%%%%%%

 [xx1 yy1 zz1]=ndgrid(origin(1)+dx(1)/2:dx(1):maximum(1)-dx(1)/2,...
                   origin(2)+dx(2)/2:dx(2):maximum(2)-dx(2)/2,...
                   origin(3)+dx(3)/2:dx(3):maximum(3)-dx(3)/2);

%figure(10)
%scatter3(origin(1)+dx(1)*(cpts(:,1)-1/2),origin(2)+dx(2)*(cpts(:,2)-1/2),origin(3)+dx(3)*(cpts(:,3)-1/2),[],cpts(:,4))
%hold on
%surf(xcoords2,ycoords2,fac1bot)
%surf(xcoords2,ycoords2,fac2bot)
%surf(xcoords2,ycoords2,fac3bot)
%surf(xcoords2,ycoords2,fac4bot)
%%surf(xcoords2,ycoords2,fac5bot)
%P=[2 1 3];
%slice(permute(xx1,P),permute(yy1,P),permute(zz1,P),permute(ID,P),dx(1)/2,dx(2)/2,1.5)
%hold off

for n=1:nreals 
 for f_loop=1:facies.n

 % statistical parameters
 alpha = 2; %(2 for normal, 0<alpha<2 for stable)
 H = 0.5;   A = 1;  % Don't know what A does
 a1 = 0.9;    a2 = .9;   a3=1.2;
 C1 = 20;    C2 = 20;   C3=1; 
 rho = 2;   
 
 xx=xx1; yy=yy1; zz=zz1;
 
 % find min and max portion of grid needed to model each facies
xx(ID~=f_loop)=[]; 
yy(ID~=f_loop)=[]; 
zz(ID~=f_loop)=[]; 
facies_min=[min(min(min(xx))) min(min(min(yy))) min(min(min(zz)))];
facies_max=[max(max(max(xx))) max(max(max(yy))) max(max(max(zz)))];
%facies_min=[1 1 min(min(min(zz)))];
%facies_max=[size(ID,1) size(ID,2) max(max(max(zz)))];
clear xx yy zz
f_min_n=floor((facies_min-origin)./dx);
f_max_n=ceil((facies_max-origin)./dx);
fm=f_max_n-f_min_n+1;
origin_facies(f_loop,:)=f_min_n;

% throw out conditioning points outside facies
fac_cpts=cpts;
for i=size(fac_cpts,1):-1:1
    try
        if (ID(fac_cpts(i,1),fac_cpts(i,2),fac_cpts(i,3))~=f_loop)
            fac_cpts(i,:)=[];
        end
    catch err
        fac_cpts(i,:)=[]; % Try/Catch in case error from cond. pt. outside domain
        cpts(i,:)=[];
        disp('conditioning point outside domain!')
    end
end

%%%% get mean and variance for each facies from conditioning points inside
%%%% that facies.
mu= mean(fac_cpts(:,4));  
means(f_loop)=mu;
%sigma = std(fac_cpts(:,4)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find stdev of each column of HRK data within each facies and average
wellcount=1;
wellstart=1;
sig=[];
for i=2:size(fac_cpts,1)
    if fac_cpts(i,1)~=fac_cpts(i-1,1) || fac_cpts(i,2)~=fac_cpts(i-1,2)
        sig(wellcount)=std(fac_cpts(wellstart:i-1,4));
        wellcount=wellcount+1;
        wellstart=i;
    end
end
sig(wellcount)=std(fac_cpts(wellstart:i,4));
sigma=mean(sig);
sigmas(f_loop)=sigma;

 % Make the Fourier kernel phi 
 Q=a1+a2+a3;
        
% k1 = repmat(-m:m,2*m+1,1);
% memory is in short supply, so don't build 3 m by m by m matrices
% should be fast enough to loop through the whole phi matrix below
 k1=-fm(1):fm(1); k2=-fm(2):fm(2); k3=-fm(3):fm(3); k1=k1*(A/fm(1)/dx(1)); 
 k2=k2*(A/fm(2)/dx(2)); k3=k3*(A/fm(3)/dx(3)); phi=ones(2*fm(1)+1,2*fm(2)+1,2*fm(3)+1);
 for i=1:2*fm(1)+1 
     for j=1:2*fm(2)+1 
         for k=1:2*fm(3)+1
 
    phi(i,j,k) = ( C1*(abs( k1(i)*theta1(1)+k2(j)*theta1(2)+k3(k)*theta1(3) ))^(rho/a1) +...
                   C2*(abs( k1(i)*theta2(1)+k2(j)*theta2(2)+k3(k)*theta2(3) ))^(rho/a2) +...
                   C3*(abs( k1(i)*theta3(1)+k2(j)*theta3(2)+k3(k)*theta3(3) ))^(rho/a3))^(1/rho); 

         end
     end
 end
        phi = phi.^(-H-Q/alpha);
         %clear k1 k2 k3;

%for k=1:2*m+1;
%  this does not matter much ...
        phi(fm(1)+1,fm(2)+1,fm(3)+1)=2*(phi(fm(1),fm(2),fm(3)+1) + phi(fm(1),fm(2)+2,fm(3)+1) + ...
        phi(fm(1)+2,fm(2)+2,fm(3)+1) + phi(fm(1)+2,fm(2),fm(3)+1) + phi(fm(1)+1,fm(2)+1,fm(3)) + ...
        phi(fm(1)+1,fm(2)+1,fm(3)))/6;
%        phi(m+1,m+1)=.1;
%end
        phi=ifftshift(phi);

        cov=real(fftn(phi.*conj(phi)))/fm(1)/fm(2)/fm(3);
        cov=fftshift(cov);    % need the covariance function for trend surf.
        
 % interpolate the trend surface of the conditioning points
 %if (docondition ==1); trendc= trend3d(cpts,cov,m,nuse); 
 if (docondition ==1)
     trendc=KrigsFast3d(fac_cpts-repmat([f_min_n 1],size(fac_cpts,1),1)+1,cov,fm);
 %%figure(1),imagesc([0,m],[0,m],trendc); hold on; plot(cpts(:,2),cpts(:,1),'o'); hold off;
 %   figure(f_loop)
 % h = slice(trendc,size(trendc,1)/2,size(trendc,2)/2,1,'noline');
 %  hold on
 %  title(strcat('Mean Field From Kriging Interpolation of Conditioning Points, facies #',num2str(f_loop)))
 %  colorbar
 %  set(h,'EdgeColor','none')
 %  scatter3(fac_cpts(:,1)-f_min_n(1)+1,fac_cpts(:,2)-f_min_n(2)+1,fac_cpts(:,3)-f_min_n(3)+1,[],fac_cpts(:,4))
 %  hold off
   
%figure(f_loop+10)
%scatter3(origin(1)+dx(1)*(fac_cpts(:,1)-1/2),origin(2)+dx(2)*(fac_cpts(:,2)-1/2),origin(3)+dx(3)*(fac_cpts(:,3)-1/2),[],fac_cpts(:,4))
%hold on
%surf(xcoords2,ycoords2,fac1bot)
%surf(xcoords2,ycoords2,fac2bot)
%surf(xcoords2,ycoords2,fac3bot)
%surf(xcoords2,ycoords2,fac4bot)
%%surf(xcoords2,ycoords2,fac5bot)
%hold off

 end

 
 % generate nreals unconditional fields
 % simulate FFT of noise field using IID normal (stable if alpha<2)
 
%for n=1:nreals 
    n=1;
 Z = randn(2*fm(1)+1,2*fm(2)+1,2*fm(3)+1);    Z = Z* ((A/fm(1)/fm(2)/fm(3))^(2/alpha));
 Z1 = randn(2*fm(1)+1,2*fm(2)+1,2*fm(3)+1);   Z1 = Z1* ((A/fm(1)/fm(2)/fm(3))^(2/alpha));
 Z = complex(Z,Z1);    % DB does not understand next step (zero mean?)
 %Z(m,m) = 0;   Z(m+1,m)=0;   Z(m,m+1)=0;    Z(m+1,m+1)=0;  
        %Z(m,m,m+1)=0;   Z(m+1,m,m+1)=0; Z(m,m+1,m+1)=0;  Z(m+1,m+1,m+1)=0;

         K = phi.*Z; % FFT of log K field obtained by multiplying FFT of noise field Z with FFT filter phi
       
        clear Z Z1;

        K = real(fftn(K))/fm(1)/fm(2)/fm(3); % Invert the FFT to get the simulated field
%        Y(1:2:end)=-Y(1:2:end); % Correct for shifting FFT integral from [-A,A] to [0,2*A] by (-1)^{i+j} 

        K=K(floor(fm(1)/2:3*fm(1)/2-1),floor(fm(2)/2:3*fm(2)/2-1),floor(fm(3)/2:3*fm(3)/2-1));    %take middle 1/4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NMonnig:  changed to mean(mean(std(K,0,3))) to get average
        % standard deviation of columns in z-direction (to match
        % variability of HRK data)
        s=mean(mean(std(K,0,3))); % average stdev for the simulated field columns
        K=(sigma/s)*(K-mean(mean(mean(K))))+mu; % convert to mean mu, standard deviation sigma
 
        if(docondition==1)  
          newpts=fac_cpts-repmat([f_min_n 1],size(fac_cpts,1),1)+1;   % make a copy of the conditioning points list so I can use the location data
          for k=1:size(fac_cpts,1)
                 % replace the random values with the new unconditioned field
          newpts(k,4)=K(newpts(k,1),newpts(k,2),newpts(k,3));             
          end
          K=K+trendc;
          %clear trendc
          %get trend surface of uncond. field
          %trendnew=trend3d(newpts,cov,m,nuse);
          trendnew=KrigsFast3d(newpts,cov,fm);
                 %make the new conditioned field
         % K=trendc+(K-trendnew);
         K=K-trendnew;
        end
    % save out facies sections of K-field
    K_facies{f_loop}=K;
    disp(strcat({'done with facies '},num2str(f_loop)))
 end
 
%%%%%%%%%% ASSEMBLE K-FIELD FROM EACH FACIES HERE %%%%%%%%%%%%%
K=zeros(m(1),m(2),m(3)); 
for i=1:m(1)
    for j=1:m(2)
        for k=1:m(3)
            facnum=ID(i,j,k);
            K(i,j,k)=K_facies{facnum}(i-origin_facies(facnum,1)+1,j-origin_facies(facnum,2)+1,k-origin_facies(facnum,3)+1);
        end
    end
end
      
%%%%%%%%%%%  plot layer by layer   %%%%%%%%%%%%%%%%%%%%
   AA=K(:,:,1);
 %  for k=1:m(3)   
 %     AA=K(:,:,k);
 % %    figure(2), imagesc(AA),pause(.1);
 % hold on
 %     figure(f_loop+1), imagesc([0,m(1)],[0,m(2)],AA'),pause(0.1)
 %     title('Wow!  Trippy Man!')
 %     colorbar
 % hold off
 %  end
   
%%%%%%%% plot 3-d view of full K-field %%%%%%%%%%%%%%%%%%%%% 
%    figure(f_loop+2)
%    g = slice(K,m(2)/2,m(1)/2,1);    
%   hold on
%   title('K Field')
%   colorbar
%   set(g,'EdgeColor','none')
%   hold off
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Look at a slice plot through 4 wells in a row
line_cpts=cpts;
for i=size(cpts,1):-1:1
    if origin(1)+dx(1)*cpts(i,1)<3 
        line_cpts(i,:)=[];
    end
end

%%%%%%%%%% Plot facies slice through 4 HRK profiles %%%%%%%%%%%%%%%%%%%
%ysurf=(origin(2)+dx(2)/2:dx(2):maximum(2)-dx(2)/2)';
%ysurf=repmat(ysurf,1,m(3));
%zsurf=(origin(3)+dx(3)/2:dx(3):maximum(3)-dx(3)/2);
%zsurf=repmat(zsurf,m(2),1);
%IDslice=reshape(ID(round((4.3-origin(1))/dx(1)),:,:),m(2),m(3));
%figure(100)
%contourf(ysurf,zsurf,IDslice)
%hold on
%scatter(origin(2)+dx(2)*(line_cpts(:,2)-1/2),origin(3)+dx(3)*(line_cpts(:,3)-1/2),[],line_cpts(:,4))
%hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
% check that conditioning points are correct
tol=10^(-10);
for i=1:size(cpts,1)
    if abs(K(cpts(i,1),cpts(i,2),cpts(i,3))-cpts(i,4))>tol
        disp('BROKEN!')
        return
    end
end
           
        fname=strcat('k3dout',int2str(n));

        % Format K values for Parflow
        K = PerFK(K);
        K = ln(K);
        % CHECK UNITS ON K!!!!! (NOTE:  PerFK takes ln(K) and returns K)
        fname=strcat('kout',int2str(n));
        dlmwrite (fname, m);
        dlmwrite (fname, K, '-append' );
        
   disp('Dunzo!')
           toc;        
end

return