%Calculates the scaling point W for the symmetric cones
function [w,w_h,l]=eval_scaling_points(problem,xc,s)
    
    %Only the scaling point for the SOCP and SDP cones
    size_w = sum(problem.soc_cones);
    w      = zeros(size_w,1);
    w_h    = zeros(sum(problem.soc_cones),1);
    l      = zeros(size_w,1);
 
    %Now calculate the scaling points and the square root of
    %the scaling point
    ie = 1;
    if(problem.n_soc_cones>0)
        %Iterate over every cone
        for(k=1:problem.n_soc_cones)
            xJx = xc(ie)^2-norm(xc(ie+1:ie+problem.soc_cones(k)-1))^2;            
            sJs = s(ie)^2 -norm( s(ie+1:ie+problem.soc_cones(k)-1))^2;            
            sx  = xc(ie:ie+problem.soc_cones(k)-1)'*s(ie:ie+problem.soc_cones(k)-1);
            w_gamma = sqrt((1+(sx/sqrt(sJs*xJx)))/2);
            %Rest of the cone
            w(ie)= sqrt(sqrt((xJx/sJs)))*(1/(2*w_gamma))*...
                                            (xc(ie)/sqrt(xJx)+1/sqrt(sJs)*s(ie)); 
           
            w(ie+1:ie+problem.soc_cones(k)-1)= sqrt(sqrt((xJx/sJs)))*(1/(2*w_gamma))*...
                                            (xc(ie+1:ie+problem.soc_cones(k)-1)/sqrt(xJx)-1/sqrt(sJs)*s(ie+1:ie+problem.soc_cones(k)-1)); 
            swJw    = sqrt(w(ie)^2-norm(w(ie+1:ie+problem.soc_cones(k)-1))^2);
            w_h(ie) = w(ie) + swJw;
            
            w_h(ie+1:ie+1+problem.soc_cones(k)-1) = w(ie+1:ie+1+problem.soc_cones(k)-1);
            w_h(ie:ie+1+problem.soc_cones(k)-1)   = sqrt((1/(2*(w(ie)+swJw))))*w_h(ie:ie+1+problem.soc_cones(k)-1);
            ie  = ie+problem.soc_cones(k);
        end
    end
end
