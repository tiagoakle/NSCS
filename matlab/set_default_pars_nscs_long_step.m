%generates the parameters structure for nscs long step with the default parameters
function pars = set_default_pars_nscs_long_step()
    pars = struct;
    %set up the default parameters
    pars.max_iter   = 200;  %maximum outer iterations
    pars.max_affine_backtrack_iter = 300;    %maximum affine backtracking steps
    pars.backtrack_affine_constant = 0.95;   %affine backtracking constant

    %xxx: changed from 0.98 for gp testing
    pars.eta           = 0.98;                %multiple of step to the boundary
    pars.stop_primal   = 1e-7;                 %stopping criteria p_res/rel_p_res<stop_primal.
    pars.stop_dual     = 1e-7;
    pars.stop_gap_res  = 1e-7;
    pars.stop_gap   = 1e-8;
    pars.stop_mu    = 1e-8;
    pars.stop_tau_kappa = 1.e-6;
    pars.solve_second_order = true;

    pars.print      = 1;                     %level of verbosity from 0 to 11
    %regularization for the linear solver
    pars.delta      = 5e-10;
    pars.gamma      = 5e-10;
    pars.max_iter_ref_rounds = 5;
    pars.neigh      = 1.0;
end

