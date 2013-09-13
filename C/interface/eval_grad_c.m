%Calls the method to evaluate the gradient
%void eval_hess_no_structs(csi k_count, csi* nK, int* tK, csi n, double* x, int* HI, int* HJ, double* HV);
function [g] = eval_grad_c(K,x)
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

    px  = full(x(permute));
    
    %Calculate the number of non zeros of the hessian
    n = size(x,1); 
    %Define the vectors for the return
    gt = libpointer('doublePtr',zeros(n,1));

    ret = calllib('nscs','eval_grad_no_structs',k_count,nK,tK,n,px,gt);
    g = zeros(n,1);
    g(permute)  = gt.value; 
    
    %void eval_hess_no_structs(csi k_count, csi* nK, int* tK, csi n, double* x, int* HI, int* HJ, double* HV);
end
