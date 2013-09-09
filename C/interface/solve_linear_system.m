function [d,CF] = solve_linear_system(H,mu,A,b,c,tau,kappa,r1,r2,r3,r4,r5,pars)
    %Calls the C code that solves the linear system

    %Debug flag 
    print_t = false;
    print_c = true;
    %Get the size of the problem
    [m,n] = size(A);

 %Set up the pointers to call the solver
 %int solve_kkt_system_no_structs(int m, int n,\
 %                               double mu,\
 %                               int* Hi,\
 %                               int* Hj,\
 %                               double* Hv,\
 %                               int Hnnz,\
 %                               int* Ai,\
 %                               int* Aj,\
 %                               double* Av,\
 %                               int Annz,\
 %                               double* b,\
 %                               double* c,\
 %                               double tau,\
 %                               double kappa,\
 %                               double delta,\
 %                               double gamma,\
 %                               double* r1,\
 %                               double* r2,\
 %                               double r3,\
 %                               double* r4,\
 %                               double r5,\
 %                               double* dy,\
 %                               double* dx,\
 %                               double* dt,\
 %                               double* ds,\
 %                               double* dk);
  
    
    [HI,HJ,HV] = find(H);
    Hnnz       = size(HI,1);
    HI(1:end)  = HI(1:end)-1;
    HJ(1:end)  = HJ(1:end)-1;
    
    [AI,AJ,AV] = find(A);
    Annz       = size(AI,1);
    AI(1:end)  = AI(1:end)-1;
    AJ(1:end)  = AJ(1:end)-1;
    
  
%    %Something is still wrong with vector c
%    % even when using a pointer
%    dc = zeros(size(c));
%    for i=1:size(c,1)
%        dc(i) = c(i);
%    end
%    c = dc;
%
%    %YUP still failing
%    db = zeros(size(b));
%    for i=1:size(b,1)
%        db(i) = b(i);
%    end
%    b = db;
%    
%    %Do r1
%    dr1 = zeros(size(r1));
%    for i=1:size(r1,1)
%        dr1(i) = r1(i);
%    end
%    r1 = dr1;
%
%
%    %Do r2
%    dr2 = zeros(size(r2));
%    for i=1:size(r2,1)
%        dr2(i) = r2(i);
%    end
%    r2 = dr2;
%
%
%
%    %Do r5
%    dr5 = zeros(size(r5));
%    for i=1:size(r5,1)
%        dr5(i) = r5(i);
%    end
%    r5 = dr5;
%
%
%    dHV = zeros(size(HV));
%    for i=1:size(HV,1)
%        dHV(i) = HV(i);
%    end
%    HV = dHV;
%
%    dAV = zeros(size(AV));
%    for i=1:size(AV,1)
%        dAV(i) = AV(i);
%    end
%    AV = dAV;
%   
%    t = HV(1);
%    HV(1) = t;
%    
%    t = AV(1);
%    AV(1) = t; 
%    
%    t = c(1);
%    c(1) = t;
%    
%    t = b(1);
%    b(1) = t;
%    
%    t = r1(1);
%    r1(1) = t;
%   
%    t = r2(1);
%    r2(1) = t;
%    
%    t = r5(1);
%    r5(1) = t;
    c          = full(c);
    b          = full(b);
    r1         = full(r1);
    r2         = full(r2);
    r5         = full(r5);

    p_c        = libpointer('doublePtr',c);
    p_b        = libpointer('doublePtr',b);
    p_r1       = libpointer('doublePtr',r1);
    p_r2       = libpointer('doublePtr',r2);
    p_r5       = libpointer('doublePtr',r5);
    p_dx       = libpointer('doublePtr',zeros(n,1));
    p_dy       = libpointer('doublePtr',zeros(m,1));
    p_ds       = libpointer('doublePtr',zeros(n,1));
    p_dt       = libpointer('doublePtr',0);
    p_dk       = libpointer('doublePtr',0);
  

    delta = 1.e-10;
    gamma = 1.e-10;
   
    %Make the call to the solver and reverse r4 and r5
    ret     = calllib('liblinear_solvers','solve_kkt_system_no_structs',m,n,mu,...
                                                                        HI,HJ,HV,Hnnz,...
                                                                        AI,AJ,AV,Annz,...
                                                                        p_b,p_c,tau,kappa,delta,gamma,...
                                                                        p_r1,p_r2,r3,p_r5,r4,...
                                                                        p_dy,p_dx,p_dt,p_ds,p_dk);
      
    d    = {}; 
    d{1} = p_dx.value;
    d{2} = p_dt.value;
    d{3} = p_dy.value;
    d{4} = p_ds.value;
    d{5} = p_dk.value;
    d{6} = 0;  %This is some kind of error flag??

    CF = [];

 end


