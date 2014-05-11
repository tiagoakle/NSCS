%Generates the parameters structure for nscs long step with the default parameters
pars = struct;
%Set up the default parameters
    pars.max_iter   = 100;  %Maximum outer iterations
    pars.max_affine_backtrack_iter = 300;    %Maximum affine backtracking steps
    pars.backtrack_affine_constant = 0.95;   %Affine backtracking constant

    %XXX: changed from 0.98 for gp testing
    pars.eta        = 0.98;                %Multiple of step to the boundary
    pars.stop_primal= 1e-5;                 %Stopping criteria p_res/rel_p_res<stop_primal.
    pars.stop_dual  = 1e-5;
    pars.stop_gap   = 1e-5;
    pars.stop_mu    = 1e-7;
    pars.stop_tau_kappa = 1.e-5;
    pars.solve_second_order = true;

    pars.print      = 1;                     %Level of verbosity from 0 to 11
    %Regularization for the linear solver
    pars.delta      = 5e-10;
    pars.gamma      = 5e-10;
    pars.max_iter_ref_rounds = 20;


