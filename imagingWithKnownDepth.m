%% To image tehe scene assuming we know the depth map

Setup_2D

%% generate the image
numDepths = 3; % this changes for each scene 
load('GT_texture2.mat') % contains reflectivity values and the range profile


depthVec = reshape(depthMap,Np*Nt,1);
% ind_set = find(depthVec~=1);
% numPoints = length(ind_set);
numPoints = Nt*Np;
ind_set = 1:Nt*Np;

Zmat = Z_test;

Z = reshape(Z_t, Np*Nt, 1); % Scene vector such that there is only one reflector in each direction
Z = Z(ind_set);
%% Build the operator

sensors_factor = [1, 0.2]; % iterates for different number of codes 

numIters = length(sensors_factor);
errors = zeros(1,numIters);

for outerLoop = 1:numIters
        
Kt = round(2*Ay*fmax/c)+1;
Kp = round(2*Ax*fmax/c)+1;


M = Kt*Kp;

Xu = zeros(numPoints,round(M*sensors_factor(outerLoop))*Kr);


if outerLoop ==1   
Phi = eye(M);
else
Phi = randn(M,round(M*sensors_factor(outerLoop)));
end

for ll = 1:Kr
 Xuf = zeros(numPoints, M);
   f = fsrc(ll);
   lambda = c/f;
   
   fay = linspace(-Ay/lambda, Ay/lambda, Kt);
   for kk = 1:Kt
   vtheta = exp(1j*2*pi*tt'*fay(kk));    
       for jj = 1:Kp
       fax = linspace(-Ax/lambda,Ax/lambda,Kp);
       
       vec = reshape(exp(1j*2*pi*tp'*fax(jj))*vtheta.', Nt*Np, 1);
       
       Xuf(:,(kk-1)*Kp+jj) = vec(ind_set).*exp(1i*2*pi*r(depthVec(ind_set))'/lambda);
       end      
   end
   Xu(:,(ll-1)*round(M*sensors_factor(outerLoop))+1:ll*round(M*sensors_factor(outerLoop))) = Xuf*Phi; 
   clear Xuf
end

Y = Xu'*Z;

x_est = real(pinvTikh(Xu',0.001)*Y);
vec_rec = zeros(Np*Nt,1);
vec_rec(ind_set) = x_est;

xMat = reshape(vec_rec,Np,Nt);

if outerLoop == 1
    x_LS = xMat;
    errors(outerLoop) = 0;
else
    errors(outerLoop) = norm(reshape(x_LS - xMat, Np*Nt,1))^2/norm(reshape(x_LS, Np*Nt,1))^2;
end
% str = strcat('img',num2str(outerLoop),'.mat');
% save(str,'xMat')


end

% save('imagingWithKnownDepth.mat','errors')
