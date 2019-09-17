%% To test CoSaMP on 2D arrays

% Xu: (Kt x Kp x Kr) x (Nt x Np x Nd)

% Y = Xu*Z; initialize x = 0;

% Compute Proxy residual: R = Y - Xu*x

% Check for stopping criteria

% Apply adjoint to R, compute teh support, take union with support pof
% previous estimate of x

% Select these columns from the matrix, solve the least squares and then
% again truncate the image

% Solving the least squares part can be done using inexpensive gradient
% steps using the proxy gradient. 

% Kp = 40, Kt = 40, Kr = 10, Nt = 78, Np = 78, Nr = 10 results in a 14.5GB
% matrix

Setup_2D

%% generate the image
numDepths = 3;
% [Zmat,rProfile,Z_t] = generateMap_2D(r,tt,tp,numDepths);
load('GT_texture2.mat')
Zmat = Z_test;
Z = reshape(Zmat, Np*Nt*Nr, 1); % Scene vector such that there is only one reflector in each direction

%% Build the operator

sensors_factor = [1 0.5 0.2 0.1 0.05];
numIters = length(sensors_factor);
errors = zeros(1,numIters);

for outerLoop = 1:numIters
    
Kt = round(2*Ay*fmax/c)+1;
Kp = round(2*Ax*fmax/c)+1;

% Pulse frequencies

Xu = zeros(Nt*Np*Nr,Kt*Kp*Kr);
M = Kt*Kp;

if outerLoop == 1
    phi = eye(M);
else
    phi = randn(M,M*sensors_factor(outerLoop));
end

for ll = 1:Kr
 Xuf = zeros(Nt*Np*Nr, M);
   f = fsrc(ll);
   lambda = c/f;
   
   fay = linspace(-Ay/lambda, Ay/lambda, Kt);
   for kk = 1:Kt
   vtheta = exp(1j*2*pi*tt'*fay(kk));    
       for jj = 1:Kp
       fax = linspace(-Ax/lambda,Ax/lambda,Kp);
       
       vec = reshape(exp(1j*2*pi*tp'*fax(jj))*vtheta.', Nt*Np, 1);
       Xuf(:,(kk-1)*Kp+jj) = reshape(vec*exp(1i*2*pi*r/lambda),Nt*Np*Nr,1);
       end      
   end
   Xu(:,(ll-1)*M*sensors_factor(outerLoop)+1:ll*M*sensors_factor(outerLoop)) = Xuf*phi; 
   clear Xuf
end

Y = Xu'*Z;

%% Modified CoSaMP for 2D arrays

x0 = zeros(size(Z));
iter = 0;
tol = 1e-3;
p = 10000*ones(size(Z));

thresh = 0.001; % threshold for being considered zero
x_inter = x0;
res = [];



while norm(p)<1e-3 || iter < 12000
    iter = iter+1;
    norm(res)
    
    % compute residual    
    res = Y - Xu'*x_inter;
    
    % compute proxy gradient 
    p = real(Xu*res/Nt/Np); % testing for real images
    if norm(p)<1e-3
        break
    end
    % identify the support candidate set
    if iter == 1
        s_x = [];
    else 
        s_x = supp(x_inter, thresh);
    end
    
    % support candidate vector
    scv = sort(unique([supp(trunc_2D(p,Np*Nt,Np*Nt))' s_x']));

    % step size obtained by solving line search analytically
    
    p_restricted = zeros(size(p));
    p_restricted(scv) = p(scv);
    
    tau = real((x_inter'*Xu - Y')*(Xu'*p_restricted)/norm(Xu'*p_restricted)^2); % enforcing tau to be real
  
    % temporary estimate
    b = zeros(size(x_inter));
    b(scv) = x_inter(scv) - tau*p(scv);
   
    
    x_inter = trunc_2D(b,Np*Nt,Np*Nt);
    
end
    
x_est = x_inter;

xMat = reshape(x_est,Np,Nt,Nr);

if outerLoop == 1
    x_LS = xMat;
    errors(outerLoop) = 0;
else
errors(outerLoop) = norm(reshape(x_LS - xMat, Np*Nt*Nr,1))^2/norm(reshape(x_LS, Np*Nt*Nr,1))^2;
end
str = strcat('img',num2str(outerLoop),'.mat');
save(str,'xMat')

end
save('CoSaMP2d_errors.mat','errors')
