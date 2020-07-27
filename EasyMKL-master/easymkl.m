%% EasMKL demo on iris dataset
% created 11-04-2018
% Okba Bekhelifi, <okba.bekhelif@univ-usto.dz>
% tic;
 function easymkl_model= easymkl(X,Y)
% load iris_dataset
% X = (LAV);
%  Y=functionaldata1.Label;

% [Y, ~] = find(irisTargets);
% X1 = X;
% Y1 = Y;
% ind = Y==1 | Y==2;
% X = X(ind, :);
% Y = Y(ind);
%  Y(Y~=1)=-1;

proptrain = 0.2;
[n, m] = size(X);
ntrain = round(n*proptrain);
ntest = n - ntrain;
rp = randperm(n);
trainset = rp(1:ntrain);
testset = rp(ntrain+1:end);

xtrain = X(trainset, :);
ytrain = Y(trainset);
xtest = X(testset, :);
ytest = Y(testset);
%% EasyMKL
gamma = 0.1;
lambda = 0.1;
rbf = @(X,Y) exp(-gamma .* pdist2(X,Y, 'euclidean').^2); % RBF kernel
d = 29; %  number of fatures in a kernel
r = 3; % number of weak kernels
[n,m] = size(xtrain);
[n1, m1] = size(xtest);
features = randi(d, [r d]);
Ks_tr = zeros(n,n,r);
for i=1:r
    tmp = xtrain(:,features(i,:));
    Ks_tr(:,:,i) = rbf(tmp, tmp);
end
Ks_ts = zeros(n1,n, r);
for i=1:r
    Ks_ts(:,:,i) = rbf(xtest(:,features(i,:)), xtrain(:,features(i,:)));
end
tracenorm = 0;
%  ytrain=functionaldata1.Label(1:48);
easymkl_model = easymkl_train(Ks_tr, ytrain', lambda, tracenorm);
 end
% acc_tr = (sum(tr_pred==ytrain)/length(ytrain))*100;
% acc_ts = (sum(ts_pred==ytest)/length(ytest))*100
% toc
% end