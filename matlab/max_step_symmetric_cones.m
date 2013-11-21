%Calculates the maximum step along the direction state.dxc, state.dtau, state.ds, state.dkappa
%to the boundary of the symmetric cones.
function a_affine = max_step_symmetric_cones(problem,state)
    %Maximum step w.r.t. the positive orthant variables
    ratios = [1;-state.xc./state.dxc;-state.s./state.ds;-state.tau/state.dtau;-state.kappa/state.dkappa];
    a_affine = min(ratios(find(ratios>1)));
    %Maximum step w.r.t the socp_cones

end    

