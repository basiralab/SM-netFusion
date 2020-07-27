clc
clear all
close all
rng default
addpath('src')
addpath('snnf')
addpath('EasyMKL-master')
addpath('libsvm-3.23/matlab')
addpath('circularGraph')
addpath('aeolianine')

%% Setting graph data simulation parameters
%%

mu1 = 0.9; % Mean value of the first Gaussian distribution
sigma1 = 0.4; % Standard deviation value of the first Gaussian distribution

mu2 = 0.7; % Mean value of the second Gaussian distribution
sigma2 = 0.6; % Standard deviation value of the second Gaussian distribution

Nf = 5; % Number of selected features

fprintf('The number of selected features is automatically set to %d.', Nf)
fprintf('\nTo change it, please set Nf variable inside run_demo.m to a different integer. \n\n');


%% Change 'displayResults' option for visualizing the learned network atlases and top selected features
%%
displayResults = 0; % input 1 if you want to visualize the estimated atlases and selected features at each run of the cross-validation algorithm
fprintf('The option for displaying the estimated atlases and selected features at each run of the cross-validation algorithm is set to %d.', displayResults)
fprintf('\nTo turn it off, please set displayResults variable inside run_demo.m to 0. \n\n');
fprintf('\nNote that displaying the results at each run will slow down the demo. \n\n');

%% Simulate graph data for running the demo
%%

% In this exemple, each class has its own statistical distribution

[data] = simulateData(mu1, sigma1, mu2, sigma2); % data samples drawn from two Gaussian distributions, each corresponding to one class

data_class1 = data.Featurematrix((data.Labels ==1),:); % retrieve samples in class 1
data_class2 = data.Featurematrix((data.Labels ==-1),:);  % retrieve samples in class 2
%%  Initialisation
k_fold=5;


addpath('C:\Users\DELL\Downloads\EasyMKL-master\EasyMKL-master')
c = cvpartition(size(data.Labels,1),'KFold',k_fold);
decision_score = zeros(size(data.Labels,1),1); %vector of decision values (indep.of treshold)
[sz1,sz2,sz3] = size(data.X);
dataFeatures = data.Featurematrix;
accuracy = zeros(size(data.Labels,1),1) ;%store accuracy from each test
predicted_Labels=[]; %vector to collect test results (Labels assigned by SV
test_Labels_vector=[]; % store ground truth in test order
ind_Nf =[];  % Store the indices of tX top discriminative features

for m = 1 : c.NumTestSets
    
    mm = num2str(m)
    
    % Create training and testing sets
    testIndex = c.test(m);
    trainIndex = c.training(m);
    train_Labels = data.Labels(trainIndex);
    train_data = data.X(trainIndex,:,:);
    test_Labels = data.Labels(testIndex);
    test_data =data.X(testIndex,:,:);
   
    
    %% SM_netFusion execution
    [Atlas1,Atlas2,topFeaturesind] = SM_netFusion(train_data,train_Labels,Nf,displayResults);
    %% Extract top Nf discriminative features
    t_C1=0;
    t_C2=0;
    for h=1:length(test_Labels)
        if (test_Labels(h)==1)
            t_C1=t_C1+1;
            test_C1(t_C1,:,:)=test_data(h,:,:);
        else
            t_C2=t_C2+1;
            test_C2(t_C2,:,:)=test_data(h,:,:);
        end
        
    end
    
    %% Extract top Nf discriminative features
    
    train_set = zeros(length(train_Labels),length(dataFeatures));
    train_Nf = zeros(length(train_Labels),Nf);
    
    % Extract the top Nf discriminative training featuress
    
    for r = 1: (length(train_Labels))
        train_subject = squeeze(train_data(r,:,:));
        train = triu(train_subject);
        train_vect = [];
        
        for i = 1: sz3
            
            for j = (i+1): sz3
                train_vect = [train_vect,train(i,j)]; % Vectorize the upper triangular part of the train matrix
            end
            
        end
        
        train_set(r,:) = train_vect; % Matrix stacking all training subjects
        
        for h = 1: Nf
            l = topFeaturesind(h);
            train_Nf(r,h) = train_set(r,l); % Discriminative training matrix
        end
        
    end
    % Extract the same ranked features from the testing network
    
    
    for r = 1: (length(test_Labels))
        test_subject = squeeze(test_data(r,:,:));
        test = triu(test_subject);
        test_vect = [];
        
        for i = 1: sz3
            
            for j = (i+1): sz3
                test_vect = [test_vect,test(i,j)]; % Vectorize the upper triangular part of the train matrix
            end
            
        end
        
        test_set(r,:) = test_vect; % Matrix stacking all testing subjects
        
        for h = 1: Nf
            l = topFeaturesind(h);
            test_Nf(r,h) = test_set(r,l); % Discriminative training matrix
        end
        
    end
    
    
    
    %% Step 3: Disease Classification using SVM classifier
    
    model = svmtrain(train_Labels,train_Nf); % Training the classifier using the training data
   [predict_Labels, accuracy, decision_values] = svmpredict(test_Labels,test_Nf,model); % Testing the classfier on the left out data (hidden/test data)
    predicted_Labels = [predicted_Labels; predict_Labels];
    test_Labels_vector = [test_Labels_vector; test_Labels];
    ind_Nf = [ind_Nf; topFeaturesind];
end


CM = confusionmat(test_Labels_vector,predicted_Labels); % Returns the confusion matrix CM determined by the known and predicted groups, respectively

True_Negative = CM(1,1);
True_Positive = CM(2,2);
False_Negative = CM(2,1);
False_Positive = CM(1,2);
Accuracy = (True_Positive + True_Negative)/(size(data.Labels,1)) * 100;
Sensitivity = (True_Positive)/(True_Positive + False_Negative) * 100;
Specificity = (True_Negative)/(True_Negative + False_Positive) * 100;
%% Display the circular graph and top features for all subjects

%size_ind_Nf = size(ind_Nf,1);
%[J,H] = scoreAcrossAllCVRuns(data,ind_Nf,Nf); %Display the circular graph of the top discriminative features
%pause(2)

%% Display final results

fprintf('\n')
disp( '                             Final results using 5-fold-CV                            ');
fprintf('\n')
disp(['****************** Accuracy = ' num2str(Accuracy) '% ******************']);
fprintf('\n')
disp(['****************** Sensitivity = ' num2str(Sensitivity) '% ******************']);
fprintf('\n')
disp(['****************** Specificity = ' num2str(Specificity) '% ******************']);

