function [AA,bb,cc,num_ter,num_var,num_con,A,constraints] = read_gp(file_name) 
    %Reads a GP in the format provided by Erling Andersen's .eo files.
    %Sets up the constraints in the form solved by coneopt:
    % 
    % The variables are ob,x_p,x_n,u,v,w.
    % [0,A,-A,-I, ,  ] [ob]    [-log(c_coef)]
    % [1,0,0 ,0 ,ix_0] [xp]    [0]
    % [0,0,0 ,0 ,IX  ] [xn]    [1]
    % [0,0,0 ,0 , I  ] [u ]    [1]
    %                  [v ]
    %                  [w ]
    %
    % num_ter: is the total number of monomial terms in the cosnstraints and objective.
    % Each "term" gives rise to a triplet u,v,w which is constrained to be in an exponential cone
    % num_var: is the number of variables in the problem
    % num_con: is the number of posinomial constraints in the problem
    %
    % AA,bb,cc define the exponential problem to be solve
    %
    % The matrix AA is constructed as follows
    % A is of size num_ter by num_var
    % IX is of size num_con by num_ter, each row corresponds to 
    % a constraint and each column to a term, the corresponding entry 
    % is 1 if the term belongs to the constraint.
    % ix_0 contains 1s in the indices that correspond to the terms in the objective.
    
    f         = fopen(file_name,'r');
    %Read the number of constraints
    num_con   = fscanf(f,'%i',1);
    %Read the number of variables
    num_var   = fscanf(f,'%i',1);
    %Read number of terms 
    num_ter   = fscanf(f,'%i',1);
    %Read the 'c' coefficients for the terms
    [c_coef,c_coef_count] = fscanf(f,'%g',num_ter);
    %Do some checking
    if(c_coef_count ~= num_ter); error('Unable to read c coefficients'); end
    %Read the constraint indices
    [constraint_index,constraint_index_count] = fscanf(f,'%i',num_ter);
    if(constraint_index_count~=num_ter); error('Unable to read the constraint indices'); end
    %Not an efficient implementation
    I = [];
    J = [];
    V = []; 
    
    %Define the vectors that indicate which terms belong to each constraint
    constraints = sparse(constraint_index+1,[1:num_ter]',ones(num_ter,1),num_con+1,num_ter,num_ter);
    
    %Define a vector for the objective
    a_0 = zeros(3*num_ter+num_var,1);
    while true
        row = fscanf(f,'%i %i %e',3);
        if isempty(row)
            break;
        end
        %The data is 0 indexed
        I = [I;row(1)+1];
        J = [J;row(2)+1];
        V = [V;row(3)];
    end
    neg_1 = sparse(num_con+1,1);
    neg_1(1) = -1;
    neg_0 = sparse(num_ter,1); %Column of zeros
    A   = sparse(I,J,V,num_ter,num_var);
    
    AA  = [[neg_0,A,-A,-speye(num_ter),sparse(num_ter,2*num_ter)];
          [neg_1,sparse(num_con+1,2*num_var+num_ter),constraints,sparse(num_con+1,num_ter)];
          [sparse(num_ter,2*num_var+2*num_ter+1),speye(num_ter)]];
    
    bb  = zeros(2*num_ter+num_con+1,1);
    bb(1:num_ter) = -log(c_coef);
    bb(num_ter+1) = 0;
    bb(num_ter+2:2*num_ter+num_con+1) = 1;
    
    cc = zeros(3*num_ter+2*num_var+1,1);
    cc(1) = 1;
    
    fclose(f);
end

