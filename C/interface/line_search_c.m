%Method to call the atd linesearch implementation
function [a,nbisections] = line_search_c(v,d,K,pars)
   
    %int linesearch_atd_no_structs( int m, int n, double*x, double*y, double*s, double tau, double kappa, double* dx,\
    %                                                                                                    double* dy,\
    %                                                                                                    double* ds,\
    %                                                                                                    double  dtau,\
    %                                                                                                    double  dkappa,\
    %                                                                                                    double  lscaff,\
    %                                                                                                    double  eta,\
    %                                                                                                    double  theta,\
    %                                                                                                    int max_backtrack,\
    %                                                                                                    int k_count,\
    %                                                                                                    int* nK,\
    %                                                                                                    int* tK,\
    %                                                                                                    double nu,\
    %                                                                                                    csi nnzH,\
    %                                                                                                    double* a,\
    %                                                                                                    int* nbacktrack);

    p_a = libpointer('doublePtr',0);
    p_i = libpointer('int32Ptr',0);

    tau = v.tau;
    kappa = v.kappa;
    
    nnzH = nnz(v.F{3}); %count the number of non zeros in H
    %Build the cone information
    
    %coneopt has a different apporach for organizing the elements of the exponential cone.
    %It concatenates all the x1s then all x2s and then all x3s
    %nscs is organized differently, for n exponential constraints we have 3n variables organized as x11,x21,x31,....,x1n,x2n,x3n

    %for now only support positive orthant and exponential cone
    % build cone:
   % K.npos = 0;
   % K.npow = 0;
   % K.nexp = N;
   % K.nlog = 0;
   % K      = getbarrpar(K);

    %The first cone is the positive orthant
    nK(1) = K.npos;
    tK(1) = 0;
    %Then the exponential cones
    nK(2:K.nexp+1) = 3;
    tK(2:K.nexp+1) = 3; 
     
    ne = 3*K.nexp; %total number of variables in the exponential cones
    %We need to permute the variables in the range (K.npos+1, to K.npos+3*K.nexp) in vectors x,s,dx,ds
    permute = zeros(ne,1);
    permute(1:3:3*K.nexp) = [1:K.nexp];
    permute(2:3:3*K.nexp) = K.nexp+[1:K.nexp];
    permute(3:3:3*K.nexp) = 2*K.nexp+[1:K.nexp]; 
    permute = [[1:K.npos];permute+K.npos];

    
    %Set the number of cones
    poscone = 0;
    if K.npos > 0; poscone = 1;end
    k_count = K.nexp+poscone; 
    m = size(v.y,1);
    n = size(v.x,1);

    px  = full(v.x(permute));
    ps  = full(v.s(permute));
    pdx = full(v.dx(permute));
    pds = full(v.ds(permute));
    y   = full(v.y);
    ret = calllib('nscs','linesearch_atd_no_structs',m,n,px,v.y,ps,v.tau,v.kappa,pdx,v.dy,pds,v.dtau,v.dkappa,pars.lscaff,pars.eta,pars.theta,pars.lsmaxit,k_count,nK,tK,K.nu,nnzH,p_a,p_i);
    a = p_a.value;
    nbisections = p_i.value;
end
