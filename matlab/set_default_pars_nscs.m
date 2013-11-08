%Set up the default parameters for nscs
    pars.max_iter   = 100;  %Maximum outer iterations
    pars.max_affine_backtrack_iter = 40;    %Maximum affine backtracking steps
    pars.max_centering_backtrack_iter = 40; %Maximum centering backtracking steps
    pars.max_c_iter = 50;                    %Maximum centering iterations per affine  iteration
    pars.backtrack_affine_constant = 0.94;   %Affine backtracking constant
    pars.backtrack_centering_constant = 0.5; %Centering backtracking constant
    pars.beta       = 0.2;                   %Stop centering when ||dx||_H<beta
    pars.theta      = 0.8;                   %Take an affine step if ||Psi^+||<theta*mu+
    %XXX: changed from 0.98 for gp testing
    pars.eta        = 0.5;                  %Multiple of step to the boundary
    pars.use_nesterov_todd_scaling = false; %Use centering points for symmetric cones
    pars.stop_primal= 1e-5;                 %Stopping criteria p_res/rel_p_res<stop_primal.
    pars.stop_dual  = 1e-5;
    pars.stop_gap   = 1e-5;
    pars.stop_mu    = 1e-5;
    pars.stop_tau_kappa = 1.e-7;
    pars.solve_second_order = true;

    pars.print      = 1;                     %Level of verbosity from 0 to 11
    %Regularization for the linear solver
    pars.delta      = 5e-10;
    pars.gamma      = 5e-10;
    pars.max_iter_ref_rounds = 2;
    pars.linear_solver = 'mixed';
    pars.centrality_measure = 1;


