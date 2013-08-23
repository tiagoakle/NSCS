#include <stdio.h>
#include <stdlib.h>
#include "smatvec.h"
#include "test_util.h"
#include "test_smatvec.h"

int tn = 5 ;
csi tAp [ ] = {0, 2, 5, 9, 10, 12} ;
csi tAi [ ] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
double tAx [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
double tb [ ] = {8., 45., -3., 3., 19.} ;
double tbt[ ] = {8.,20.,13.,6.,17} ;
double tx [ ] = {1., 2. ,  3., 4.,  5.} ;

void LoadSystemData(char * nameA,\
                    char * namen,\
                    char * namem,\
                    csi ** Ap,\
                    csi ** Ai,\
                    double** Av,\
                    csi *m,\
                    csi *n,\
                    csi *nnz,\
                    double** nvector,\
                    double** mvector)
{
   //XXX: If the files are made in int and we read size_t we will have a crash
    read_csv_size(nameA,m,n,nnz);
    //printf("File %s contains sparse matrix of size %i, %i ,%i\n",nameA,*m,*n,*nnz);
    csi* I = calloc(*nnz,sizeof(csi));
    csi* J = calloc(*nnz,sizeof(csi));
    double* V = calloc(*nnz,sizeof(double));
    read_csv_triplets(nameA,I,J,V);
    
    //printf("First non zero %d, %d ,%lf\n",*I,*J,*V);
    *nvector = calloc(*n,sizeof(double));
    read_vector(namen,*nvector);
    *mvector = calloc(*m,sizeof(double));
    read_vector(namem,*mvector);

    (*Ap) = calloc(*n+1,sizeof(csi));
    (*Ai) = calloc(*nnz,sizeof(csi));
    (*Av) = calloc(*nnz,sizeof(double));
    
    //After reading the test matrix convert it into CSR

    int status;    
    status  =  umfpack_di_triplet_to_col(*m,*n,*nnz,I,J,V,*Ap,*Ai,*Av,NULL);
    //printf("Pack status %i\n",status); 
    //printf("First non zero packed %d, %d, %lf\n",**Ap,**Ai,**Av);
    free(I);
    free(J);
    free(V);
 
}

START_TEST (small_example)
{
  
  csi m = 5;
  csi n = 5;
  csi nnz = 12;
   
  double* y = calloc(m,sizeof(double));
  //Now call the smvp
  double alpha = 1;
  dspmv(m, n, alpha, tAp, tAi, tAx, y, tx);
  //Calculate the difference between the stored result and the 
  //calculated
  
  double t = 0;
  double acc = 0;
  int i = 0;
  for(i=0;i<m;i++)
  {
      t = (y[i]-tb[i]);
      acc += t*t;
  }

    printf("Dist sq %lf\n",acc);
    ck_assert(acc<1.e-5);
    //At the end free 
    free(y);
}
END_TEST


START_TEST (small_example_transpose)
{
  
  csi m = 5;
  csi n = 5;
  csi nnz = 12;
   
  double* y = calloc(m,sizeof(double));
  //Now call the smvp
  double alpha = 1;
  dsTpmv(m, n, alpha, tAp, tAi, tAx, y, tx);
  //Calculate the difference between the stored result and the 
  //calculated
  
  double t = 0;
  double acc = 0;
  int i = 0;
  for(i=0;i<m;i++)
  {
      t = (y[i]-tbt[i]);
      acc += t*t;
      //  if(i<5) printf("y[%i]-tb[%i]: %lf\n",i,i,t);
  }

    printf("Dist sq %lf\n",acc);
    //At the end free 
    free(y);
    ck_assert(acc<1.e-5);
}
END_TEST


START_TEST (direct_product_no_empty_cols)
{
  csi n,m,nnz;
  csi* Ai, * Ap;
  double* Av; 
  double*x, *yc;

  char* nameA = "./test/test_data/Matrix1.csv";
  char* namex = "./test/test_data/x1.csv";
  char* namey = "./test/test_data/y1.csv";
  LoadSystemData(nameA,namex,namey,&Ap,&Ai,&Av,&m,&n,&nnz,&x,&yc);

  double* y = calloc(m,sizeof(double));
  //Now call the smvp
  double alpha = 1;
  dspmv(m, n, alpha, Ap, Ai, Av, y, x);
  //Calculate the difference between the stored result and the 
  //calculated
  
  double t = 0;
  double acc = 0;
  int i = 0;
  for(i=0;i<m;i++)
  {
      t = (y[i]-yc[i]);
      acc += t*t;
  }
    write_vector_to_csv("./test/test_data/result1.csv",y,m);
    printf("Dist sq %lf\n",acc);
    ck_assert(acc<1.e-5);
    //At the end free 
    free(Ap);
    free(Ai);
    free(Av);
    free(y);
    free(yc);
    free(x);
}
END_TEST

START_TEST (direct_product_with_empty_cols)
{
  csi n,m,nnz;
  csi* Ai, * Ap;
  double* Av; 
  double*x, *yc;

  char* nameA = "./test/test_data/Matrix3.csv";
  char* namex = "./test/test_data/x3.csv";
  char* namey = "./test/test_data/y3.csv";
  LoadSystemData(nameA,namex,namey,&Ap,&Ai,&Av,&m,&n,&nnz,&x,&yc);

  double* y = calloc(m,sizeof(double));
  //Now call the smvp
  double alpha = 1;
  dspmv(m, n, alpha, Ap, Ai, Av, y, x);
  //Calculate the difference between the stored result and the 
  //calculated
  
  double t = 0;
  double acc = 0;
  int i = 0;
  for(i=0;i<m;i++)
  {
      t = (y[i]-yc[i]);
      acc += t*t;
  }
    write_vector_to_csv("./test/test_data/result1.csv",y,m);
    printf("Dist sq %lf\n",acc);
    ck_assert(acc<1.e-5);
    //At the end free 
    free(Ap);
    free(Ai);
    free(Av);
    free(y);
    free(yc);
    free(x);
}
END_TEST


START_TEST (transpose_product_no_empty_cols)
{

    //Here we test A'x = y
    
  csi n,m,nnz;
  csi* Ai, * Ap;
  double* Av; 
  double*nvec, *mvec;

  char* nameA = "./test/test_data/Matrix1.csv";
  char* namen = "./test/test_data/y2.csv";
  char* namem = "./test/test_data/x2.csv";
     
  LoadSystemData(nameA,namen,namem,&Ap,&Ai,&Av,&m,&n,&nnz,&nvec,&mvec);


    double* prod_res = calloc(n,sizeof(double));
    //Now call the smvp
    double alpha = 1;
    dsTpmv( m, n, alpha, Ap, Ai, Av, prod_res, mvec);
    //Calculate the difference between the stored result and the 
    //calculated
    
    double t = 0;
    double acc = 0;
    int i = 0;
    for(i=0;i<n;i++)
    {
        t = (prod_res[i]-nvec[i]);
        acc += t*t;
    }
    printf("Dist sq %lf\n",acc);
    ck_assert(acc<1.e-5);
    //At the end free 
    free(Ap);
    free(Ai);
    free(Av);
    free(prod_res);
    free(mvec);
    free(nvec);
}
END_TEST

START_TEST (transpose_product_with_empty_cols)
{

    //Here we test A'x = y
    
  csi n,m,nnz;
  csi* Ai, * Ap;
  double* Av; 
  double*nvec, *mvec;

  char* nameA = "./test/test_data/Matrix3.csv";
  char* namen = "./test/test_data/y4.csv";
  char* namem = "./test/test_data/x4.csv";
     
  LoadSystemData(nameA,namen,namem,&Ap,&Ai,&Av,&m,&n,&nnz,&nvec,&mvec);


    double* prod_res = calloc(n,sizeof(double));
    //Now call the smvp
    double alpha = 1;
    dsTpmv( m, n, alpha, Ap, Ai, Av, prod_res, mvec);
    //Calculate the difference between the stored result and the 
    //calculated
    
    double t = 0;
    double acc = 0;
    int i = 0;
    for(i=0;i<n;i++)
    {
        t = (prod_res[i]-nvec[i]);
        acc += t*t;
    }
    printf("Dist sq %lf\n",acc);
    ck_assert(acc<1.e-5);
    //At the end free 
    free(Ap);
    free(Ai);
    free(Av);
    free(prod_res);
    free(mvec);
    free(nvec);
}
END_TEST


START_TEST (scalar_case)
{

    //Here we test A'x = y with A = 1, x = 2 
    
  csi n,m,nnz;
  csi* Ai, * Ap;
  double* Av; 
  double*nvec, *mvec;

  char* nameA = "./test/test_data/Matrix0.csv";
  char* namen = "./test/test_data/y0.csv";
  char* namem = "./test/test_data/x0.csv";
     
  LoadSystemData(nameA,namen,namem,&Ap,&Ai,&Av,&m,&n,&nnz,&nvec,&mvec);


    double* y = calloc(n,sizeof(double));
    //Now call the smvp
    double alpha = 1;

    //printf("A: %lf, x: %lf, y: %lf \n",*Av,*nvec,*mvec);

    dspmv(m, n, alpha, Ap, Ai, Av, y, mvec);
    //Calculate the difference between the stored result and the 
    //calculated
    
    double t = 0;
    double acc = 0;
    int i = 0;
    for(i=0;i<n;i++)
    {
        t = (y[i]-nvec[i]);
        acc += t*t;
    }
    printf("Dist sq %lf\n",acc);
    ck_assert(acc<1.e-5);
    //At the end free 
    free(Ap);
    free(Ai);
    free(Av);
    free(y);
    free(mvec);
    free(nvec);
}
END_TEST

Suite* matvec_suite(void)
{
    Suite* suite = suite_create("Smatvec");
    TCase *tc = tcase_create("Case1");
    tcase_add_test(tc,small_example);
    tcase_add_test(tc,small_example_transpose);
    tcase_add_test(tc,scalar_case);
    tcase_add_test(tc,direct_product_no_empty_cols);
    tcase_add_test(tc,direct_product_with_empty_cols);
    tcase_add_test(tc,transpose_product_no_empty_cols);
    tcase_add_test(tc,transpose_product_with_empty_cols);
    suite_add_tcase(suite,tc);
    return suite;
}

