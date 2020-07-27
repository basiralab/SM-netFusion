%% Main function of SM_netFusion framework for a fast and accurate classification.
% Details can be found in the original paper:
% Islem Mhiri and Islem Rekik. "Supervised Multi-topology Network
% Cross-diffusion for Population-driven Brain Network Atlas Estimation"
%


%   ---------------------------------------------------------------------

%     This file contains the implementation of three key steps of our SM_netFusion framework:
%     (1) Class-specific feature extraction and clustering,
%     (2) Class-specific supervised multi-topology network cross-diffusion  and  
%     (3) Discriminative connectional biomarker identification:
%
%                 [AC1,AC2,ind] = SM_netFusion(train_data,train_Labels,Nf,displays)
%
%                 Inputs:
%
%                          train_data: ((n/5) × 4) × m × m) tensor stacking the symmetric matrices of the training subjects
%                                      n the total number of subjects
%                                      m the number of nodes
%
%                          train_Labels: ((n/5) × 4) × 1) vector of training labels (e.g., -1, 1)
%
%                          Nf: Number of selected features
%
%                          displays: Boolean variables [0, 1].
%                                    if displays = 1 ==> display(Atlas of group 1, Atlas of group 2, top features matrix and the circular graph)
%                                    if displays = 0 ==> no display
%                 Outputs:
%                         AC1: (m × m) matrix stacking the atlas of group 1
%
%                         AC2: (m × m) matrix stacking the atlas of group 2
%
%                         ind: (Nf × 1) vector stacking the indices of the top disciminative features
%
%
%     To evaluate our framework we used Leave-One-Out cross validation strategy.



%To test SM-netFusion on random data, we defined the function 'simulateData' where the size of the dataset is chosen by the user.
% ---------------------------------------------------------------------
%     Copyright 2019 Islem Mhiri, Sousse University.
%     Please cite the above paper if you use this code.
%     All rights reserved.
%     """

%%------------------------------------------------------------------------------






function [AC1,AC2,ind] = SM_netFusion(train_data,train_Labels,Nf,displays)
s1 = 0;
s2 = 0;
[sz1,sz2,sz3] = size(train_data);
for h = 1: length(train_Labels)
    
    if (train_Labels(h) == 1)
        s1 = s1+1;
        XC1(s1,:,:) = squeeze(train_data(h,:,:));
        
    else
        s2 = s2+1;
        XC2(s2,:,:) = squeeze(train_data(h,:,:));
    end
    
end


%Disentangling the heterogeneous distribution of the input_ networks using SIMLR clustering method

% C1 group
k = [];

for l = 1: s1
    k1 = squeeze((XC1(l,:,:)));
    k2 = k1(:); %vectorize the matrix
    k = [k;k2'];
end

[t, S1, F1, ydata1,alpha1] = SIMLR(k,2,2);

%C2 group
kk = [];

for l = 1: s2
    kk1 = squeeze(abs(XC2(l,:,:)));
    kk2 = kk1(:); %vectorize the matrix
    kk = [kk;kk2'];
end
[t1, S2, F2, ydata2,alpha2] = SIMLR(kk,2,2);


% After using SIMLR, we extract each cluster independently for both classes

% C1 group
qC11 = 1;
qC12 = 1;


for qC1 = 1: s1
    
    if t(qC1) == 1
        Ca1(qC11,:,:) = abs(XC1(qC1,:,:));
        La1(qC11)=1;
        qC11 = qC11+1;
        
    elseif t(qC1) == 2
        Ca2(qC12,:,:) = abs(XC1(qC1,:,:));
        La2(qC12)=-1;
        qC12 = qC12+1;
        
        
    end
    
end

%C2 group
qC21 = 1;
qC22 = 1;


for qC2 = 1: s2
    
    if t1(qC2) == 1
        Cn1(qC21,:,:) = (XC2(qC2,:,:));
        Ln1(qC21)=1;
        qC21 = qC21+1;
        
    elseif t1(qC2) == 2
        Cn2(qC22,:,:) = (XC2(qC2,:,:));
        Ln2(qC22)=-1;
        qC22 = qC22+1;
        
        
    end
    
end

% For each cluster, we non-linearly diffuse and fuse all networks into a local centered network atlas using SNFCall

% Setting all the parameters.
K = 20;%number of neighbors, usually (10~30)
alpha = 0.5; %hyperparameter, usually (0.3~0.8)
T = 20; %Number of Iterations, usually (10~20)


for l = 1: (qC11-1)
    
    ll = num2str(l);
    Datap1.(['datap',ll,'']) = Standard_Normalization(squeeze(Ca1(l,:,:)));
    Distp1.(['distp',ll,'']) = dist2( Datap1.(['datap',ll,'']), Datap1.(['datap',ll,'']));
    Wp1.(['Wp1',ll,'']) = affinityMatrix(Distp1.(['distp',ll,'']), K, alpha);
end

Wall1 = struct2cell(Wp1);


for l = 1: (qC12-1)
    ll = num2str(l);
    Datap2.(['datap',ll,'']) = Standard_Normalization(squeeze(Ca2(l,:,:)));
    Distp2.(['distp',ll,'']) = dist2( Datap2.(['datap',ll,'']), Datap2.(['datap',ll,'']));
    Wp2.(['Wp1',ll,'']) = affinityMatrix(Distp2.(['distp',ll,'']), K, alpha);
end

Wall2 = struct2cell(Wp2);
Walla=[Wall1;Wall2];
La=[La1,La2];

% C2 group
for l = 1: (qC21-1)
    ll = num2str(l);
    Datan1.(['datan',ll,'']) = Standard_Normalization(squeeze(Cn1(l,:,:)));
    Distn1.(['distn',ll,'']) = dist2( Datan1.(['datan',ll,'']), Datan1.(['datan',ll,'']));
    Wn1.(['Wn1',ll,'']) = affinityMatrix(Distn1.(['distn',ll,'']), K, alpha);
end

Walln1 = struct2cell(Wn1);


for l = 1: (qC22-1)
    ll = num2str(l);
    Datan2.(['datan',ll,'']) = Standard_Normalization(squeeze(Cn2(l,:,:)));
    Distn2.(['distn',ll,'']) = dist2( Datan2.(['datan',ll,'']), Datan2.(['datan',ll,'']));
    Wn2.(['Wn1',ll,'']) = affinityMatrix(Distn2.(['distn',ll,'']), K, alpha);
end

Walln2 = struct2cell(Wn2);

Walln=[Walln1;Walln2];
Ln=[Ln1,Ln2];

%% SNF_all

AC1 = SNF_all(Walla,K,La,T); % Global network atlas for C1 group

AC2 =SNF_all(Walln,K,Ln,T);% Global network atlas for C2 group
%% Step 1: Estimation of a centered and representative input_ network atlas
% Extract C1 group and C2 group


D = abs(AC1-AC2);% Difference between matrices

%% Step 2:  Discriminative connectional biomarker identification % % % % %
D = triu(D,1); %Upper triangular part of matrix
D1 = D(find(D));
D1 = D1.';
D2 = sort(D1,'descend');% Ranking features
Dif = D2(1:Nf); % Extract the top Nf ranked features
ind = [];

for i = 1: Nf
    [a,b,pos] = intersect(Dif(i),D1); % Select the indices of top Nf ranked features
    ind = [ind;pos];
end
qq = [];
qq1 = [];

for h = 1 : Nf
    X = ind(h);
    [q,q1] = map_index_to_position_in_matrix(X,sz3);
    qq = [qq;q];
    qq1 = [qq1;q1];
end

topFeatures = zeros(sz3);

for i = 1: Nf
    topFeatures(qq(i),qq1(i)) = D2(i);
    topFeatures(qq1(i),qq(i)) = D2(i);
end

%% Display Atlas 1, Atlas 2 and the circular graph of the top discriminative features


if (displays)
    
    
    imagesc(AC1),title('Atlas of group 1 ','Color','b') % Display Atlas 1 of each LOO-CV iteration
    
    pause(2)
    
    imagesc(AC2) ,title('Atlas of group 2 ','Color','b') % Display Atlas 2 of each LOO-CV iteration.
    
    pause(2)
    
    imagesc(topFeatures) ,title('Top discriminative features ','Color','b') % Display Atlas 2 of each LOO-CV iteration.
    
    pause(2)
    
    close()
    
    J = circularGraph(topFeatures)
    
    pause(2)
    
end
end





