function wire_vector_to_csv(name, v)
    n = size(v,1);
    fid = fopen(name,'w');
    fprintf(fid,'%d\n',n);
    for i=1:n
        fprintf(fid,'%f\n',full(v(i)));
    end
    fclose(fid);
end

