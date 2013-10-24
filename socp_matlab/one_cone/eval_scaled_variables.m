%Calculate the scaling point for the pair x,s
%Defined as the point here H(w)x = s
%Evaluate H(w)
%Calculate W = H^{1/2}(w) and Wi=H^{-1/2}(w)
%lambda = W*x = Wi*s

%parameters:
%n, size of the cone,
%x, primal variable
%s, dual variable
function [g,H,W,Wi,lambda] = eval_scaled_variables(n,x,s)
    J   = diag(sparse([1;-ones(n-1,1)]));
   
    %Scaling points
    Jx = J*x;
    xJx = x'*Jx;
    Js = J*s;
    sJs = s'*Js;
    sx  = s'*x;
    w_gamma = sqrt((1+(sx/sqrt(sJs*xJx)))/2);
    w      = sqrt(sqrt((xJx/sJs)))*(1/(2*w_gamma))*(x/sqrt(xJx)+1/sqrt(sJs)*J*s); 
    w2     = sqrt(sqrt((sJs/xJx)))*(1/(2*w_gamma))*(J*x/sqrt(xJx)+1/sqrt(sJs)*s);
    
    %Gradient at x
    g      = -1/(xJx)*Jx;
 
   %Evaluate the hessian 
    Jw  = J*w;
    wJw = w'*Jw;
    H   = 2/(wJw)^2*(Jw*Jw') - 1/(wJw)*J; 

    %Evaluate u = w^{1/2}
    u    = zeros(n,1);
    swJw = sqrt(wJw);
    u(1) = w(1) + swJw;
    u(2:n) = w(2:n);
    u    = u*sqrt((1/(2*(w(1)+swJw))));

    
    %Evaluate W
    Ju  = J*u;
    uJu = u'*Ju;
    W   = 2/(uJu)^2*(Ju*Ju') - 1/(uJu)*J; 
 
    %Evaluate u^-1
    u = 1/uJu*Ju;
    Ju = J*u;
    uJu = u'*Ju;

    %Evaluate W^-1 
    Ju   = J*u;
    uJu  = u'*Ju;
    Wi   = 2/(uJu)^2*(Ju*Ju') - 1/(uJu)*J; 

    lambda = W*x;
 

end

