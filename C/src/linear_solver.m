%This function calls the dynamically loadable
%library to solve a linear system

function x = linear_solver(A,b)
    %We assume that the library with name liblinear_solvers has been loaded
    [I,J,V] = find(A); %Convert to triplet form
    n       = size(A,1);
    nnz     = n*2;
    lpX = libpointer('doublePtr',zeros(n,1));
    ret     = calllib('liblinear_solvers','solve_linear_system',lpX,I,J,V,nnz,b,n);
    x       = lpX.Value;
end

