%Writes a matlab matrix into a binary file 
function write_matrix_bin(name, v) 
    f = fopen(name,'w');
    [I,J,V] = find(v); %transform to coo format
    [m,n] = size(v);
    nnz = size(I,1);
    fwrite(f,m,'int32');
    fwrite(f,n,'int32');
    fwrite(f,nnz,'int32');
    fwrite(f,I,'int32');
    fwrite(f,J,'int32');
    fwrite(f,V,'double');
end
