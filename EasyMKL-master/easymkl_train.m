function [model] = easymkl_train(Ks, y,lambda, tracenorm)
%EASYMKL_TRAIN train an EasyMKL model[1]
% Input : 
%         Ks :  Set of Kernels
%         y  :  [double] Training labels 1|-1
%         lambda : regularization parameter [0,1]
%         tracenorm : 0|1 logical value whether to normalize trace or not
% Output :
%        model : EasyMKL model Struct
%              .gamma   : [Nx1] [double] instance coefficients 
%                         N : number of training examples  
%              .bias    : 1x1 double bias
%              .weights : [1xL] [double] kernel weights
%              .labels  : [NxN] [double] training labes in diagonal format
% Requirements:
% -  MOSEK : quadprog function 
% References:
% [1] Fabio Aiolli and Michele Donini 
%      EasyMKL: a scalable multiple kernel learning algorithm
%      Paper @ http://www.math.unipd.it/~mdonini/publications.html
% created 11-04-2018
% last modfied -- -- --
% Okba Bekhelifi, <okba.bekhelif@univ-usto.dz>

[sz1,sz2,sz3] = size(Ks);

kk = [];
 for l = 1: sz1
     kk1 = squeeze(Ks(l,:,:));
     kk2 = kk1(:); %vectorize the matrix
     kk = [kk;kk2'];
 end
Ks_tr = multipleK(kk,sz3);
nr_kernels = size(Ks_tr, 3);
% trace normalization

if(tracenorm)
    for i=1:nr_kernels
        Ks_tr(:,:,i) = (Ks_tr(:,:,i)*size(Ks_tr(:,:,i),1)) / trace(Ks_tr(:,:,i));
    end
end
% sum of kernels
K = sum_kernels(Ks_tr);
% 
x = optimize(K, y, lambda);
YY = diag(y);
bias = 0.5 * x' * K * YY * x;
yg = x'.*y;
weights = zeros(1,nr_kernels);
for i=1:nr_kernels
    weights(i) = yg*Ks_tr(:,:,i)*yg';
end
weights = weights ./ sum(weights);
K = sum_kernels(Ks_tr, weights);
% 
x = optimize(K, y, lambda);
% model
model.gamma = x;
model.bias = bias;
model.weights = weights;
model.labels = YY;
end
function Kernels = multipleK(x,sz1)


N = size(x,1);
KK = 0;
sigma = [((sz1/2)+0.5):-0.5:1];
Diff = (dist2(x,x));
[T,INDEX]=sort(Diff,2);
[m,n]=size(Diff);
allk = 1; 
ti=1;

    if allk < (size(x,1)-1)
        TT=mean(T(:,2:(allk+1)),2)+eps;
        Sig=(repmat(TT,1,n)+repmat(TT',n,1))/2;
        Sig=Sig.*(Sig>eps)+eps;
        for j = 1:length(sigma)
            W=normpdf(Diff,0,sigma(j)*Sig);
            Kernels(:,:,KK+ti) = (W + W')/2;
            ti = ti+1;
        end
    end
end
function [x] = optimize(K, y, lambda)
YY = diag(y);
KLL = (1-lambda) .* YY .* K .* YY;
LID = diag(lambda* ones(1,length(y)));
Q = 2 * (KLL+LID);
p = zeros(length(y),1);
G = - diag(ones(length(y),1));
h = zeros(size(K,1),1);
A = double([y<0;y>0]);
tt=size(A);
b = [1;1];
t=size(b);
[x, fval, exitflag,output] = quadprog(Q,p,G,h,A,b);
end
