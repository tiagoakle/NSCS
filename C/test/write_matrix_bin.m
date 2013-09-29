%Writes a matlab matrix into a csvary file 
function write_matrix_bin(name, v) 
    f = fopen(name,'w');
    [I,J,V] = find(v); %transform to coo format
    [m,n] = size(v);
    nnz = size(I,1);
    fprintf('Saving maxtrix of size m %i n %i ,nnz %i\n',m,n,nnz);
    fwrite(f,int32(m),'int32');
    fwrite(f,int32(n),'int32');
    fwrite(f,int32(nnz),'int32');
    fwrite(f,int32(I-1),'int32');
    fwrite(f,int32(J-1),'int32');
    fwrite(f,V,'double');
    fclose(f);
end
