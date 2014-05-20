function [ret,problem]=validate_problem_structure(problem,varargin)
    %If varargin not empty then it contains a fid for the prints
    if(~isempty(varargin))
        fid = varargin{1};
    else
        %print to stdout 
        fid = 1;
    end

    ret = 0;
    %Make sure the number of constrained and free variables adds to n;
    tot_constraints = problem.n_free+problem.n_pos+...
                      problem.n_exp_cones*3;
    
    if(tot_constraints ~= problem.n)
        fprintf(fid, 'Error, sum of all constraints must equal n\n');
        ret = -8;
        return;
    end
       
    %Make sure A,b,c are defined 
    if(~isfield(problem,'A'))
        fprintf(fid, 'Constraint matrix A is not defined\n');
        ret = -8;
        return;
    end
 
    if(~isfield(problem,'b'))
        fprintf(fid,'Vector b is not defined\n');
        ret = -8;
        return;
    end

    if(~isfield(problem,'c'))
        fprintf(fid,'Vector c is not defined\n');
        ret = -8;
        return;
    end
    
    %Check the dimensions of b and c
    if(any(size(problem.c)~=[problem.n,1]))
        fprintf(fid,'objective must be a nx1 vector \n');
        ret = -8;
        return; 
    end
 
    %Check the dimensions of b and c
    if(any(size(problem.b)~=[problem.m,1]))
        fprintf(fid,'b must be a mx1 vector \n');
        ret = -8;
        return;
    end

    [m,n]=size(problem.A);
    if(m~=problem.m)
        fprintf(fid,'A must have problem.m rows\n');
        ret = -8;
        return;
    end
    
    if(n~=problem.n)
        fprintf(fid,'A must have problem.n cols\n');
        ret = -8;
        return;
    end

  
end
