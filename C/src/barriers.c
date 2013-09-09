#include "barriers.h"
#include "common.h"
/**
 * Returns true if x is primal feasible.
 * @param prob the problem definition strucutre
 * @param x the present primal point
 * @return true if x is feasible
 */
bool primal_feas(problem_t prob, double* x)
{
    bool feas = true;
    csi i = 0;
    csi var_ix = 0;
    //Iterate over all the cones and add the respective complexities 
    for(i = 0; i < prob.k_count; i++)
    {
        feas = cone_primal_feas(prob.tK[i],x+var_ix,prob.nK[i]); 
        var_ix += prob.nK[i];
        if(!feas) break;
    }
    return feas;
}

/**
 * Returns true if x is dual feasible.
 * @param prob the problem definition strucutre
 * @param x the present dual point
 * @return true if x is feasible
 */
bool dual_feas(problem_t prob, double*x )
{
    bool feas = true;
    csi i = 0;
    csi var_ix = 0;
    //Iterate over all the cones and add the respective complexities 
    for(i = 0; i < prob.k_count; i++)
    {
        feas = cone_dual_feas(prob.tK[i],x+var_ix,prob.nK[i]); 
        var_ix += prob.nK[i];
        if(!feas) break;
    }
    return feas;
}


/**
 * Calculates the complexity parameter for 
 * the barrier of the constraint cone. 
 * @param prob problem definition structure
 * @return complexity of the barrier
 */
double barrier_complexity(problem_t prob)
{
    double nu = 0.;
    csi nnzH = 0;
    csi i = 0;
    //Iterate over all the cones and add the respective complexities 
    for(i = 0; i< prob.k_count; i++)
    {
        nu += cone_barrier_complexity(prob.tK[i],prob.nK[i]); 
    }
    return nu;
}

/**
 * Evaluates the gradient of the barrier at x
 * @param prob problem definition
 * @param x    present primal point
 * @param grad vector to store the gradient
 */
void eval_grad(problem_t prob, double* x, double* grad)
{
    double nu = 0.;
    csi val_ix = 0;
    csi i = 0;
    //Iterate over all the cones and evaluate the respective gradients
    for(i = 0; i < prob.k_count; i++)
    {
        cone_barrier_grad(prob.tK[i],grad+val_ix,x,prob.nK[i]);
        val_ix += prob.nK[i];
    }
}

/**
 * Evaluates the Hessian of the barrier at x
 * @param prob problem definition
 * @param x    present primal point
 * @param H    matrix to store the hessian
 */
void eval_hess(problem_t prob, double* x, state_t state)
{
    double nu = 0.;
    csi val_ix = 0;
    csi offset = 0;
    csi i = 0;
    //Iterate over all the cones and evaluate the respective gradients
    for(i = 0; i < prob.k_count; i++)
    {
        //cone_barrier_fills HI,HJ,HV with the coo format of the hessian
        //and return the number of non zeros in the matrix block
        offset+=cone_barrier_hessian(prob.tK[i],\
                             state.H.I+offset,\
                             state.H.J+offset,\
                             state.H.V+offset,\
                             x+val_ix,\
                             prob.tK[i],\
                             prob.delta);
        val_ix += prob.nK[i];
    }
}

/**
 * Calculates the number of non zeros for the
 * Hessian of the barrier of the constraint cone.
 * @param prob problem definition structure
 */
csi hessian_nnz(problem_t prob)
{
    csi nnzH = 0;
    csi i = 0;
    //Iterate over all the cones 
    for(i = 0; i< prob.k_count; i++)
    {
        nnzH += cone_nnz(prob.tK[i],prob.nK[i]); 
    }
    return nnzH;
}

//Calls the function that returns the nnz for the cone type
csi cone_nnz(int type, csi n)
{
    switch(type)
    {
        case 0:
         return pos_orthant_nnz(n);
        break;
        case 1:
        case 2:
        case 3:
        case 4:
        break;
    }
    return 0;
}

/*
 * Returns the complexity for the barrier of a cone 
 * of type (type) and size (dim)
 * @param the type of cone 
 * @param the size of the cone
 */

int cone_barrier_complexity(int type, csi n)
{
    switch(type)
    {
        case 0:
         return pos_orthant_complexity(n);
        break;
        case 1:
        case 2:
        case 3:
        case 4:
        break;
    }
    return 0;
}

/*
 * Returns true if the x is feasible for the dual cone 
 * of the cone with type (type) and size (dim)
 * @param the type of cone 
 * @param the size of the cone
 */

bool cone_dual_feas(int type, double* x, csi n)
{
    switch(type)
    {
        case 0:
            return pos_orthant_feas( x, n); 
        break;
        case 1:
        case 2:
        case 3:
        case 4:
        break;
    }
    return 0;
}

/*
 * Returns true if the x is feasible for the primal cone 
 * of the cone with type (type) and size (dim)
 * @param the type of cone 
 * @param the size of the cone
 */

bool cone_primal_feas(int type, double* x, csi n)
{
    switch(type)
    {
        case 0:
            return pos_orthant_feas( x, n); 
        break;
        case 1:
        case 2:
        case 3:
        case 4:
        break;
    }
    return 0;
}


/*
 * Evaluates the barrier function for the given index type
 * @param type cone type
 * @param x point at which the barrier is evaluated
 * @param n the size of x
 */

double cone_barrier_val(int type, double*x, csi n)
{
    switch(type)
    {
        case 0:
         return pos_orthant_val(x,n);
        break;
        case 1:
        case 2:
        case 3:
        case 4:
        break;
    }
    return 0;

}

/*
 * Evaluates the gradient for the barrier of the cone 
 * @param type the type of the cone
 * @param g the vector where the gradient will be stored
 * @param variable_ix the index of the first variable of the cone
 * @param x point at which the hessian is evaluated
 * @param n number of variables in the cone
 */
void cone_barrier_grad(int type, double*g, double*x, csi n)
{
    switch(type)
    {
        case 0:
             pos_orthant_grad(g,x,n);
        break;
        case 1:
        case 2:
        case 3:
        case 4:
        break;
    }

}

/**
 * Evaluates the hessian of the barrier
 * @param type the type of cone
 * @param HI the vector of row indices of the non zeros
 * @param HJ the vector of col indices of the non zeros
 * @param HV the vector of non zeros
 * @param x point where the gradient is evaluated
 * @param n the number of variables in the cone
 * @param delta regularization
 * @return the number of non zeros
 */
csi cone_barrier_hessian(int type, int* HI, int* HJ, double* HV, double*x, csi n, double delta)
{
    switch(type)
    {
        case 0:
            return  pos_orthant_hessian(HI,HJ,HV,x,n,delta);
        break;
        case 1:
        case 2:
        case 3:
        case 4:
        break;

    }

}

//Returns the complexity of the barrier
csi pos_orthant_complexity(csi n)
{
    return n;
}

//Returns the number of non zeros of the hessian of the positive orthant
csi pos_orthant_nnz(int n)
{
   return n; 
}

//Evaluates the barrier of the positive orthant
double pos_orthant_val(double* x, csi n)
{    
    double f = 0;
    csi i = 0;
    for(i=0;i<n;i++)
    {
        if(x[i] < 0)
        {
            f = INFINITY;
            break;
        }
        f += -log(x[i]);
    }
    return f;
}

//Returns true if x is feasible in the positive orthant
bool pos_orthant_feas(double* x, csi n)
{    
    bool feas = true;
    csi i = 0;
    for(i=0;i<n;i++)
    {
        if(x[i]<0){feas = false; break;}
    }
    return feas;
}

//Evaluates the gradient of the positive orthant
void pos_orthant_grad(double* grad, double *x, csi n)
{
    int i = 0;
    for(i=0;i<n;i++)
    {
        grad[i] += -1./x[i];
    }

}
//Evaluates the hessian of the positive orthant
csi pos_orthant_hessian(int* HI, int* HJ, double* HV, double*x, csi n, double delta)
{ 
    csi i = 0;
    for(i=0;i<n;i++)
    {
        HI[i] = i;
        HJ[i] = i;
        HV[i] = 1./x[i]*1./x[i] + delta; 
    }

    return n;
}


