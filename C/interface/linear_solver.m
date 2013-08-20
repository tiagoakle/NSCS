%This function calls the dynamically loadable
% library to solve a linear system
function x = linear_solver(A,b)
    %We assume that the library with name liblinear_solvers has been loaded
    [I,J,V] = find(A); %Convert to triplet form
    I = I-1;           %Shift to c notation
    J = J-1;
    n       = size(A,1);
    nnz     = size(V,1); %We are using a dense matrix
    %Define the pointer for the return  
    lpX     = libpointer('doublePtr',zeros(1,n));
    ret     = calllib('liblinear_solvers','solve_linear_system',lpX,I,J,V,nnz,b,n);
    x       = lpX.Value;
end

