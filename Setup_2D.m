% Setup for fast forward and adjoint transform

c =  3e8;%299792458;
% min and max excitation frequencies
fmin = 2e9;
fmax = 4e9;

% A = aperture size (in meters), 
%     set as WLA wavelengths of average freq, in both axes
WLAx = 19;
Ax = WLAx*c/2/fmax;%WLAx*c*2/(fmin+fmax);


WLAy = 19;
Ay =  WLAy*c/2/fmax;%WLAy*c*2/(fmin+fmax);


R0 = 20;
R = 3;

Nr = 10;

r = linspace(R0,R0+R,Nr);

theta0 = sin(pi/4)/2;
phi0 = sin(pi/4)/2;

% grids for tau / r  domain
% number of tau samples --- "Nyquist" is A*fmax/c
Nt = 78;%round(4*Ax*fmax/c);
tt = linspace(-theta0, theta0, Nt);  % 0.4330 = sin(pi/3)/2;
deltt = mean(diff(tt));

Np = 78;%round(4*Ay*fmax/c);
tp = linspace(-phi0, phi0, Np);
deltp = mean(diff(tp));

Kr = 5;

fsrc = linspace(fmin, fmax, Kr);

% number of samples to take on the omega_tau_1 and \omega_tau_p axis
%   (i.e. number of array elements) 
%   (lambda_min/2 spacing would be 2*A*fmax/c)
%   (not much seems to change when increased passed critical)



% numCodes = Kt*Kp;
% code = eye(Kt*Kp);