% This function compute the score of each top features vector and display the circular graph
function [J,H] = scoreAcrossAllCVRuns(input_data,ind_selectedFeatures,Number_topFeatures)
% Initialisation

ind_Nf = ind_selectedFeatures.';
C = ind_Nf(:);
C = C';
ind_Nfuni = unique(C);
[sz1,sz2,sz3] = size(input_data.X);
[sz1_Nf,sz2_Nf] = size(ind_selectedFeatures);

for k = 1 : length(ind_Nfuni)
    kk = num2str(k);
    [i,j] = find(ind_Nf == ind_Nfuni(k));
    i = sz2_Nf-i;
    ii.(['i',kk]) = i;
    jj.(['j',kk]) = j;
end

for m = 1 : length(ind_Nfuni)
    mm = num2str(m);
    score_vector(m) = sum(ii.(['i',mm]))/sz1_Nf;
end
S = sort(score_vector);
SS = S((length(ind_Nfuni)-(Number_topFeatures-1)):length(ind_Nfuni));
sv = [];
for i = 1 : Number_topFeatures
    [a,b,pos] = intersect(SS(i),score_vector);
    sv = [sv;pos];
end
y =[];
qq = [];
qq1 = [];
for i = 1 : Number_topFeatures
    y = [y;ind_Nfuni(sv(i))];
    
end

for h = 1 : Number_topFeatures
    X = y(h);
    [q,q1] = map_index_to_position_in_matrix(X,sz3);
    qq = [qq;q];
    qq1 = [qq1;q1];
end
H = zeros(sz3);
for i = 1 : Number_topFeatures
    H(qq(i),qq1(i)) = SS(i);
    H(qq1(i),qq(i)) = SS(i);
    
end
J = circularGraph(H)
end