SetupCommon % script to set up all the parameters

%% Target Image
img = imread('Shapes.png'); % file containing ground truth reflectivity points
img = imresize(rgb2gray(img),[Np,Nt]);
img = double(img);
Z = double(reshape(img,Nt*Np,1));

Z = reshape(repmat(Z,1,Nr),Np*Nt*Nr,1) ;


%% Build the operator


sensors_factor = [0.01]; % to determine number of codes. contains fractions to multiply M with
numIters = length(sensors_factor);
errors = zeros(1,numIters);

    
Kt = round(2*Ay*fmax/c)+1; % Number of array elements along axis 1
Kp = round(2*Ax*fmax/c)+1; % Number of array elements along axis 2

Xu = zeros(Nt*Np*Nr,Kt*Kp*Kr); % Array operator declaration
M = Kt*Kp; % total number of array elements

% Build the array operator 

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
   Xu(:,(ll-1)*M+1:ll*M) = Xuf; 
   clear Xuf
end


%% Conventional Imaging

Xu_sep = [real(Xu)'; imag(Xu)'];
Y_full = Xu_sep*Z;

% noise = randn(size(Y_full));
% noise = noise/norm(noise)*norm(Y_full)*0.01;

% Y_full = Y_full+noise;% uncomment to add noise

% delta = 5; % uncomment to set regularization parameter

% [Xu_sep_pinv,Uf,Sf,Vf] = pinvTikh(Xu_sep, delta);% uncomment to regularize

X_rec_full = pinv(Xu_sep)*Y_full;
% noise_rec_full = Xu_sep_pinv*noise;% uncomment if noise is present
%save('X_rec_full', 'X_rec_full'); % uncomment if you need to store data

%opSNR(1) = 20*log10(norm(X_rec_full)/norm(noise_rec_full)); % uncomment
%under noise

%X_rec_full_noisy = X_rec_full+noise_rec_full; % uncomment under noise
%% Aperture coding

%delta_d = delta*ones(numel(sensors_factor)); % regularization for coding 

for ii = 1:numIters
    
    numCodes = M*sensors_factor(ii);
    phi = randn(M,numCodes);
    
    Xuc = zeros(Nt*Np*Nr, numCodes*Kr); % The aperture coded operator
    
    for jj = 1:Kr
        Xuc(:,(jj-1)*numCodes+1:jj*numCodes) = Xu(:,(jj-1)*M+1:jj*M)*phi;
    end
    
    Xuc_sep = [real(Xuc)'; imag(Xuc)'];
    
    % Generate the aperture coded measurements by using the new operator
    
    Y_c = Xuc_sep*Z;
    
%     uncomment the nex five lines under noisy regimes
%     noise = randn(size(Y_c));
%     noise = noise/norm(noise)*norm(Y_c)*0.1;

%     Y_c = Y_c+noise;
    
%    [Xuc_sep_pinv,Ufc,Sfc,Vfc] = pinvTikh(Xuc_sep, 5); 

%     X_rec_c = Xuc_sep_pinv*Y_c;
    X_rec_c = pinv(Xuc_sep)*Y_c; % comment under noisy regime

    
    errors(ii) = norm(X_rec_full - X_rec_c)^2/norm(X_rec_full)^2;
    
    str = strcat('X_rec_',num2str(sensors_factor(ii)),'.mat');
%     save(str,'X_rec_c') ; % uncomment if reconstruction needs to be
%     stored
    
    
end
        
%     save('const_range.txt','errors') % uncomment if errors need to be
%     stored
    