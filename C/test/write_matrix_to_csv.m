function write_matrix_to_csv(name, A)
    %Writes a csv file name which contains the IJV representation of A
    % the first row contains the values n,m,nnz
    % the ith row contains I[i]-1,J[i]-1,V[i].
    [m,n] = size(A);
    nnzA   = nnz(A);
    [I,J,V] = find(A);

    fid = fopen(name,'w');
    fprintf(fid,'%d,%d,%d\n',m,n,nnzA);
    for i=1:n
        fprintf(fid,'%f,%f,%f\n',full(I(i)-1),full(J(i)-1),full(V(i)));
    end
    fclose(fid);
end

