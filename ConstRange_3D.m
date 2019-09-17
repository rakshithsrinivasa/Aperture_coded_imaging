SetupCommon

%% Target Image
%img = imread('Shapes.png');
%img = imresize(rgb2gray(img),[Np,Nt]);
%img = double(img);
%Z = double(reshape(img,Nt*Np,1));

load('img3d.mat');
Z = reshape(img,Nt*Np*Nr,1);
%% Build the operator


sensors_factor = [0.01];
numIters = length(sensors_factor);
errors = zeros(1,numIters);
opSNR = zeros(1,numIters+1);
    
Kt = round(2*Ay*fmax/c)+1;
Kp = round(2*Ax*fmax/c)+1;

Xu = zeros(Nt*Np*Nr,Kt*Kp*Kr);
M = Kt*Kp;


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
delta = 0.01;
Xu_sep = [real(Xu)'; imag(Xu)'];
Y_full = Xu_sep*Z;

noise1 = randn(size(Y_full));
noise1 = noise1/norm(noise1)*norm(Y_full)*0;

Y_full = Y_full;

[Xu_sep_pinv,Uf,Sf,Vf] = pinvTikh(Xu_sep,delta);
%ranks(1) = sum(diag(Sf)>1e-5);
X_rec_full = Xu_sep_pinv*Y_full;
noise_rec_full = Xu_sep_pinv*noise1;

opSNR(1) = 20*log10(norm(X_rec_full)/norm(noise_rec_full)); 

X_rec_full_noisy = X_rec_full+noise_rec_full;

save('X_rec_full_noisy_3d.mat','X_rec_full_noisy')

%% Aperture coding

delta_d = delta*ones(numel(sensors_factor));%2*[15,36,65,149];
Id = eye(M);
for ii = 1:numIters
    
        numCodes = M*sensors_factor(ii);
	rng('shuffle');
	col_ind = randperm(M);
        %phi = randn(M,nu,numCodes);
	phi = Id(:,col_ind(1:numCodes)); 
        phi(1) 
    Xuc = zeros(Nt*Np*Nr, numCodes*Kr);
    
    for jj = 1:Kr
        Xuc(:,(jj-1)*numCodes+1:jj*numCodes) = Xu(:,(jj-1)*M+1:jj*M)*phi;%/sqrt(numCodes);
    end
    
    Xuc_sep = [real(Xuc)'; imag(Xuc)'];
    Y_c = Xuc_sep*Z;
    
    %noise2 = randn(size(Y_c));
    %noise2 = noise2/norm(noise2)*norm(Y_c)*0.01;
    noise2 = zeros(size(Y_c));
    for jj = 1:Kr
	noise2((jj-1)*numCodes+1:jj*numCodes,1) = phi'*noise1((jj-1)*M+1:jj*M,1);%/sqrt(numCodes);
    end
    disp('norm of compressed noise is:')
    norm(noise2)
    disp('norm of Y_c is:')
    norm(Y_c)

    Y_c = Y_c;
    delta_hat = delta_d(ii);
    [Xuc_sep_pinv,Ufc,Sfc,Vfc] = pinvTikh(Xuc_sep,delta_hat);

    X_rec_c = Xuc_sep_pinv*Y_c;
    noise_rec_c = Xuc_sep_pinv*noise2;
    
    errors(ii) = norm(X_rec_full - X_rec_c)^2/norm(X_rec_full)^2;
    opSNR(ii+1) = 20*log10(norm(X_rec_c)/norm(noise_rec_c));

    X_rec_c_noisy = X_rec_c+noise_rec_c;    

    str = strcat('X_rec_noisy',num2str(sensors_factor(ii)),'_3d.mat');
    save(str,'X_rec_c_noisy')
   end
    save('const_range_errors_noisy_3d.mat','errors')
    %save('opSNR_constRange.mat','opSNR') 
    

    
