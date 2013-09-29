%Evaluates the build residuals function
function[p_res,d_res,g_res,n_p_res,n_d_res,n_g_res,rel_gap] = build_residuals_c(K,y,x,s,tau,kappa,A,b,c,p_relstop,d_relstop,g_relstop)
    
    %Build the permutation for A,x,s,c
    ne = 3*K.nexp; %total number of variables in the exponential cones
    %We need to permute the variables in the range (K.npos+1, to K.npos+3*K.nexp) in vectors x,s,dx,ds
    permute = zeros(ne,1);
    permute(1:3:3*K.nexp) = [1:K.nexp];
    permute(2:3:3*K.nexp) = K.nexp+[1:K.nexp];
    permute(3:3:3*K.nexp) = 2*K.nexp+[1:K.nexp]; 
    permute = [[1:K.npos];permute+K.npos]';

    px  = full(x(permute));
    ps  = full(s(permute));
    pc   = full(c(permute));
    
    %Permute A and build coo format
    pA    = A(:,permute);

  %  px  = full(x);
  %  ps  = full(s); 
  %  pc   = full(c);

    [Ai,Aj,Av] = find(pA);
    Ai = Ai-1;
    Aj = Aj-1;
    [m,n] = size(A);
    nnzA   = nnz(A);

    %Define the work vectors
    psi = libpointer('doublePtr',zeros(n,1));
    hpsi = libpointer('doublePtr',zeros(n,1));
    cent = libpointer('doublePtr',0); 

    p_p_res = libpointer('doublePtr',zeros(m,1));
    p_d_res = libpointer('doublePtr',zeros(n,1));
    p_g_res = libpointer('doublePtr',0);

    pn_p_res = libpointer('doublePtr',0);
    pn_d_res = libpointer('doublePtr',0);
    pn_g_res = libpointer('doublePtr',0);
    prel_gap = libpointer('doublePtr',0);

    calllib('nscs','calculate_residuals_no_structs',...
    Ai,Aj,Av,m,n,nnzA,...
    full(b),pc,...
    full(y),px,tau,ps,kappa,p_relstop,d_relstop,g_relstop,...
    p_p_res,p_d_res,p_g_res,...
    pn_p_res,pn_d_res,pn_g_res,prel_gap);
    
    %Assign the return values
    p_res = p_p_res.value;
    d_res = zeros(n,1);
    d_res(permute) = p_d_res.value;
    

    %'from build_residuals_c' 
    %norm(p_d_res.value-pA'*y-ps+tau*pc)
    %norm(d_res - A'*y-s+tau*c)

    g_res = p_g_res.value;
    n_p_res = pn_p_res.value;
    n_d_res = pn_d_res.value;
    n_g_res = pn_g_res.value;
    rel_gap = prel_gap.value;
     
%void calculate_residuals_no_structs(csi* Ai, csi* Aj, double* Av, csi m, csi n, csi nnz,
%                                    double* b, double* c,\
%                                    double* y, double* x, double tau, double* s, double kappa,\
%                                    double p_relstop, double d_relstop, double g_relstop,\
%                                    double* p_res, double* d_res, double* g_res,\
%                                    double* n_p_res, double* n_d_res, double* n_g_res, double* rel_gap)
end
