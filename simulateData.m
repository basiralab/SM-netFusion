%simulateData function using Gaussian distribution
function [data] = simulateData(mu1, sigma1, mu2, sigma2)
rng(1);
Featurematrix = [];
Labels = [];
X =[];
data =[];
prompt = 'Select the number of class 1 graphs: ';
C1 = input(prompt)
while (isempty(C1) == 1)
     prompt = 'Please choose a number: ';
    C1 = input(prompt)
end
while (C1<5)
    prompt = 'Please choose a number >4: ';
    C1 = input(prompt)
end 

prompt = 'Select the number of class 2 graphs: ';
C2 = input(prompt)
while (isempty(C2) == 1)
     prompt = 'Please choose a number: ';
    C2 = input(prompt)
end
while (C2<5)
    prompt = 'Please choose a number >4: ';
    C2 = input(prompt)
end 


prompt = 'Select the number of nodes (i.e., ROIS for brain graphs): ';
m = input(prompt)
while (isempty(m) == 1)
    prompt = 'Please choose a number >20: ';
    m = input(prompt)
end 
while ((m<21) == 1) 
    prompt = 'Please choose a number >20: ';
    m = input(prompt)
end 



N = C1+C2;
dataC1 = normrnd(mu1,sigma1,[C1,m,m]);% Normal random ditribution
dataC2 = normrnd(mu2,sigma2,[C2,m,m]);% Normal random ditribution
data1 = [dataC1;dataC2];
% %% Drawing samples from two different distributions to simulate both classes
% h1 = histogram(dataC1)
% hold on
% h2 = histogram(dataC2)

for i = 1:N
data1(i,:,:)=squeeze(data1(i,:,:))-diag(diag(squeeze(data1(i,:,:)))); % Eliminate self symetry (diagonal=0)
data1(i,:,:) = (squeeze(data1(i,:,:))+(squeeze(data1(i,:,:)))')./2; % Insure data symetry
 t = triu(squeeze(data1(i,:,:)),1); % Upper triangular part of matrix
 x = t(find(t)); % Vectorize the triangle
 x1 = x.';
 Featurematrix = [Featurematrix;x1];
 data.Featurematrix = Featurematrix;
 data.X = data1;
end
data.Labels = [ones(C1,1);-1*ones(C2,1)]; %  Define labels
 end 