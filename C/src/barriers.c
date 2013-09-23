#include "barriers.h"
#include "common.h"
#include <stdio.h>
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
 * Evaluates the value of the primal barrier at x
 * @param prob the problem definition
 * @param x the point at which to evaluate
 * @return the value of the barrier
 */
double eval_barrier(problem_t prob, double* x)
{
    double barrier_val = 0.;
    csi i = 0;
    csi var_ix = 0;
    for(i = 0;i<prob.k_count;i++)
    { 
        barrier_val += cone_barrier_val(prob.tK[i],x+var_ix,prob.nK[i]);
        var_ix += prob.nK[i];
    }
    return barrier_val;
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
        cone_barrier_grad(prob.tK[i],grad+val_ix,x+val_ix,prob.nK[i]);
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
    //index of the first variable in this cone
    csi var_ix = 0;
    //non zeros calculated already
    csi offset = 0;
    csi i = 0; //iterator
    csi row_shift = 0; //first index of this hessian sub_block in the hessian
    //Iterate over all the cones and evaluate the respective gradients
    for(i = 0; i < prob.k_count; i++)
    {
        //cone_barrier_fills HI,HJ,HV with the coo format of the hessian
        //and return the number of non zeros in the matrix block
        offset+=cone_barrier_hessian(prob.tK[i],\
                             state.H.I+offset,\
                             state.H.J+offset,\
                             state.H.V+offset,\
                             row_shift,\
                             x+var_ix,\
                             prob.nK[i],\
                             prob.delta);
        var_ix += prob.nK[i];
        row_shift = var_ix;
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
            return exp_nnz();
            break;
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
            return exp_complexity();
            break;
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
            return exp_dual_feas(x);
            break;
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
            return exp_primal_feas(x);
            break;
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
            return exp_val(x);
            break;
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
            exp_grad(g,x);
            break;
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
 * @param row_shift constant to add to each entry in HI and HJ
 * @param x point where the gradient is evaluated
 * @param n the number of variables in the cone
 * @param delta regularization
 * @return the number of non zeros
 */
csi cone_barrier_hessian(int type, int* HI, int* HJ, double* HV,csi row_shift, double*x, csi n, double delta)
{
    switch(type)
    {
        case 0:
            return  pos_orthant_hessian(HI,HJ,HV,row_shift,x,n,delta);
            break;
        case 1:
        case 2:
        case 3:
            return exp_hessian(HI,HJ,HV,row_shift,x,delta);
            break;
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
        grad[i] = -1./x[i];
    }

}
//Evaluates the hessian of the positive orthant
csi pos_orthant_hessian(int* HI, int* HJ, double* HV, csi row_shift, double*x, csi n, double delta)
{ 
    csi i = 0;
    for(i=0;i<n;i++)
    {
        HI[i] = i+row_shift;
        HJ[i] = i+row_shift;
        HV[i] = 1./x[i]*1./x[i] + delta; 
    }

    return n;
}

//--------------------------Exponential cone barrier functions

//Returns the complexity of the barrier
csi exp_complexity()
{   
    return 3;
}

//Returns the number of non zeros of the hessian of the positive orthant
csi exp_nnz()
{
   return 9; 
}

//Evaluates the barrier of the positive orthant
double exp_val(double* x)
{    
  double  x1    = x[0];
  double  x2    = x[1];
  double  x3    = x[2];

  double  logx2 = log(x2);
  double  logx3 = log(x3);
  return  -log(x3*(logx2-logx3)-x1)-logx2-logx3; 
}

//Returns true if x is feasible in the positive orthant
bool exp_primal_feas(double* x)
{    
    bool feas = true;
    double logx2 = log(x[1]);
    double logx3 = log(x[2]);
    double tmp1  = logx2-logx3;
    double psi   = x[2]*tmp1 - x[0];
    // Feasibility:
    if ((psi < 0) || (x[1]<0) || (x[2]<0))
        feas = false;
        
    return feas;
}

//Returns true if x is feasible in the positive orthant
bool exp_dual_feas(double* x)
{   
   
    double dual[3];
    dual[0] = -x[2];
    dual[1] = exp(1)*x[1];
    dual[2] = -x[0];
    //return (u<0)&&(-u*exp(v/u-1)<w); //This is Boyds exam criteria
    return exp_primal_feas(dual); //This is Skajaas criteria
}

//Evaluates the gradient of the positive orthant
void exp_grad(double* grad, double *x)
{
  double  x1    = x[0];
  double  x2    = x[1];
  double  x3    = x[2];

  double  logx2 = log(x2);
  double  logx3 = log(x3);
  double  x2m1  = 1./x2;
  double  x3m1  = 1./x3;
  double  tmp1  = logx2-logx3;
  double  psi   = x3*tmp1 - x1;
  double  psim1 = 1./psi;
  double  xi    = tmp1 - 1;
  
  grad[0] = psim1;
  grad[1] = -x2m1*(x3*psim1 + 1);
  grad[2] = -xi*psim1 - x3m1;

}

//Evaluates the hessian of the positive orthant
csi exp_hessian(int* HI, int* HJ, double* HV,csi row_shift, double*x, double delta)
{ 
    double  x1    = x[0];
    double  x2    = x[1];
    double  x3    = x[2];
  
    double  logx2 = log(x2);
    double  logx3 = log(x3);
    double  x2m1  = 1./x2;
    double  x3m1  = 1./x3;
    double  tmp1  = logx2-logx3;
    double  psi   = x3*tmp1 - x1;
    double  psim1 = 1./psi;
    double  xi    = tmp1 - 1;
    
    double psi2  = psi*psi;
    double psim2 = psim1*psim1;
    double x2m2  = x2m1*x2m1;
    double x3m2  = x3m1*x3m1;
  
    double el11  = psim2;
    double el21  = -x3*x2m1*psim2;
    double el31  = -xi*psim2;
    double el22  = psim2*x2m2*(x3*psi + x3*x3 + psi2);
    double el32  = psim2*x2m1*(x3*xi - psi);
    double el33  = psim2*(x3m1*psi + xi*xi + psi2*x3m2);
     
    HV[0] = el11+delta;
    HV[1] = el21;
    HV[2] = el31;
    HV[3] = el21;
    HV[4] = el22+delta;
    HV[5] = el32;
    HV[6] = el31;
    HV[7] = el32;
    HV[8] = el33+delta;
    
    HI[0]  = 0+row_shift;
    HI[1]  = 1+row_shift;
    HI[2]  = 2+row_shift;
    HI[3]  = 0+row_shift;
    HI[4]  = 1+row_shift;
    HI[5]  = 2+row_shift;
    HI[6]  = 0+row_shift;
    HI[7]  = 1+row_shift;
    HI[8]  = 2+row_shift;

    HJ[0]  = 0+row_shift;
    HJ[1]  = 0+row_shift;
    HJ[2]  = 0+row_shift;
    HJ[3]  = 1+row_shift;
    HJ[4]  = 1+row_shift;
    HJ[5]  = 1+row_shift;
    HJ[6]  = 2+row_shift;
    HJ[7]  = 2+row_shift;
    HJ[8]  = 2+row_shift;

    return 9;
}
//---------------------End of exponential cone barrier 
