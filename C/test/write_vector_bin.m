%Writes a matlab vector into a binary file 
function write_vector_bin(name, v)
    f = fopen(name,'w');
    fwrite(f,v,'double');
    fclose(f);
end
