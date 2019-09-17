% Setup for fast forward and adjoint transform

% Initialize required parameters

c =  3e8;

% min and max excitation frequencies
fmin = 2e9;
fmax = 4e9;

% Ax, Ay = aperture size (in meters), 
%     set as WLA wavelengths of the max freqeuncy
WLAx = 39;
Ax = WLAx*c/2/fmax;%WLAx*c*2/(fmin+fmax);


WLAy = 39;
Ay =  WLAy*c/2/fmax;%WLAy*c*2/(fmin+fmax);


R0 = 20; % lower range limit
R = 0.5*c/fmax; % range extent
Rmax = R0 + R; % upper range limit

Nr = 1; % number of points in along the range axis % replace with range profile if known

r = linspace(R0, Rmax, Nr); % uniform sampling along the range axis

theta0 = sin(pi/2)/2; % Angular extent along elevation 
phi0 = sin(pi/2)/2; % Angular extent along azimuthal 

% grids for tau / r  domain
% number of tau samples --- "Nyquist" is A*fmax/c
Nt = round(4*Ax*fmax/c); % Discretization in angular domain (elevation)
tt = linspace(-theta0, theta0, Nt);  
deltt = mean(diff(tt));

Np = round(4*Ay*fmax/c); % Discretization in angular domain (azimuthal)
tp = linspace(-phi0, phi0, Np);
deltp = mean(diff(tp));

Kr = 10; % number of excitation frequencies 

fsrc = linspace(fmin, fmax, Kr); % excitation frequencies
