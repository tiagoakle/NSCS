%Evaluates the centrality measure
function [centmeas] = eval_centmeas_c(K,xa,sa,mua);
    poscones = 0;
    if(K.npos > 0)
        poscones = 1;
    end
    k_count = K.nexp + poscones;
    
    cones = 0;
    %The first cone is the positive orthant
    if(K.npos ~= 0)
        nK(1) = K.npos;
        tK(1) = 0;
        cones = 1;
    end
    %Then the exponential cones
    nK(1+cones:K.nexp+cones) = 3;
    tK(1+cones:K.nexp+cones) = 3; 
     
    ne = 3*K.nexp; %total number of variables in the exponential cones
    %We need to permute the variables in the range (K.npos+1, to K.npos+3*K.nexp) in vectors x,s,dx,ds
    permute = zeros(ne,1);
    permute(1:3:3*K.nexp) = [1:K.nexp];
    permute(2:3:3*K.nexp) = K.nexp+[1:K.nexp];
    permute(3:3:3*K.nexp) = 2*K.nexp+[1:K.nexp]; 
    permute = [[1:K.npos];permute+K.npos]';

    pxa  = full(xa(permute));
    psa  = full(sa(permute));
   
    %Calculate the number of non zeros of the hessian
    nnzH = 9*K.nexp+K.npos;
    n    = size(xa,1); 
    
    %Define the work vectors
    psi = libpointer('doublePtr',zeros(n,1));
    hpsi = libpointer('doublePtr',zeros(n,1));
    cent = libpointer('doublePtr',0); 

    HI = libpointer('int32Ptr',zeros(nnzH,1));
    HJ = libpointer('int32Ptr',zeros(nnzH,1));
    HV = libpointer('doublePtr',zeros(nnzH,1));

    ret = calllib('nscs','eval_cent_meas_no_structs',k_count,nK,tK,0,pxa,psa,n,nnzH,mua,psi,hpsi,cent);
    centmeas = cent.value;
end
%double eval_cent_meas_no_structs(csi k_count, csi* nK, int* tK, double delta, double* x, double* s, csi n, csi nnzH, double* psi, double * hpsi)
