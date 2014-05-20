%Reads a sol file and extracts the objective and the primal solution
function [prim_vars, obj] = read_mosek_sol(file_path)
    f = fopen(file_path,'r');
    file_string = char(fread(f))';
    %Find the objective 
    tok=regexp(file_string,'OBJECTIVE\s*:\s*([\d|\.|e|\-|\+]*)','tokens');
    obj=tok{1}{1};
    obj=str2double(obj);
    %Where do the primal variables start
    st=regexp(file_string,'PRIMAL\s*VARIABLES\s*\nINDEX\s*ACTIVITY\s*\n');
    en=regexp(file_string,'DUAL\s*VARIABLES');
    file_string = file_string(st+42:en-7);
    vars = regexp(file_string,'(\d+\s+[\d|\.|e|\-|\+]*)','tokens');
    prim_vars = [];
    %Extract each var
    for j=1:size(vars,2)
        val = regexp(vars{j}{1},'\d+\s+([\d|\.|e|\-|\+]*)','tokens');
        prim_vars = [prim_vars;str2double(val{1}{1})];
    end
    fclose(f);
end
