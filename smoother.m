%*********************************************************%
% Universal Approach Smoother Function                    %
% Developed from D. Pedretti and D. Fernandez-Garcia 2013 %
% Coded by Hayden Fischer, Fall 2013                      %
%*********************************************************%

clear all
disp('Beginning smoother...');


%**********************************
% Import files and define variables
%**********************************

% Import times and densities
fid   = load('well.extract.txt');
t     = fid(:,2);           % X-axis values
n     = fid(:,3);           % Particle concentration
steps = size(t,1);          % Number of time steps

% Define remaining variables
C     = 0.0;  % Functional form of P for a datum
alpha = 0.9;  % Adaptive factor
h_0   = 40;   % Sampling lag interval, bandwidth (hr)
hmin  = 1.5;  % minimum allowable bandwidth (to avoid delta functions at sub-collection time intervals)


%****************
% Generate arrays
%****************

p_0t    = zeros(steps,1);
lambda  = zeros(steps,1);      % Normalized local data density
h_1t    = zeros(steps,1);
h_2t    = zeros(steps,1);
p_1t    = zeros(steps,1);      % Particle probability density (PDF)


%****************
% Solve variables
%****************

%  for simplicity of notation, I will loop over j for making the
%  calculations of densities, and loop over i for the particle arrival
%  times

for j = 1:steps
    term=10^-10;
    for i=1:steps
        term = term + (n(i)/sqrt(2*pi*h_0^2) * exp( -(t(j)-t(i))^2 / (2*h_0^2)));
    end
    p_0t(j) = term;
end

G = exp(sum(log(p_0t)) / steps);

for i = 1:steps
    lambda(i) = (p_0t(i) / G)^-alpha;
    h_1t(i) = max([hmin   h_0 * lambda(i)]);
end

for j=1:steps
    for i=1:steps
       p_1t(j) = p_1t(j) + n(i)/sqrt(2*pi*h_1t(i)^2) * exp(-(t(j)-t(i))^2 / (2*h_1t(i)^2));  
    end
end


%*************
% Plot results
%*************

figure(2)
plot(t, n, 'o', 'Color', [0.7 0.7 0.7])
hold on
plot(t, p_1t, 'k')
legend('p_1t','SLIM data')
title('Relative Concentration vs. Time')
xlabel('Extraction Time [hours]')
ylabel('Particle Concentration [ppm]')
hold off

smoothbtc = p_1t(51:460);
save('smoothbtc.txt','-ascii','smoothbtc')

exit