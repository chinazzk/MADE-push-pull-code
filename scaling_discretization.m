disp('Beginning scaling_discretization...')
tic;

% Import K field
K = importdata('Kout3.sa');
N = 128 * 128 * 640;
K = K(2:N+1,1);
K = reshape(K, [128 128 640]);

%K = K./24;  % changing from meter/day to meter/hour
header = [128 128 640];

x = 1; y = 1; z = 1;
i=3;
%for i = 0:5
    n = int2str(i);
    a = int2str(x); b = int2str(y); c = int2str(z);
    disp(strcat('Writing madeperm',n,'.sa'))
    disp(strcat('Scaling: [',a,',',b,',',c,']'))
    %madeperm = upscaledK(K, x, y, z);
    madeperm = PerFK(K);
    title = strcat('Fbm',n,'.sa');
    dlmwrite(title, header, '\t');
    dlmwrite(title, madeperm, '-append');
    toc;
   % x = x*2; y = y*2; z = z*2;  % Increase scaling values
%end

%exit

% 4   4   20    0
% 8   8   40    1
% 16  16  80    2
% 32  32  160   3
% 64  64  320   4
% 128 128 640   5