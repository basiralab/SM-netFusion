% This function map the index into a position in a matrix
function[i,j]=map_index_to_position_in_matrix(X,sz)


k = 1;
for i = 1 : sz
 for j = (i+1): sz
     Ms(i,j) = k;
     k = k+1;

 end
end

[i,j] = find(Ms == X);
end
    
    
    
    