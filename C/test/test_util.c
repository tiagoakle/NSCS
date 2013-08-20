#include "test_util.h"
int dump_matrix_to_csv(char* name, int* I, int* J, double* V, int n)
{
    int i=0;
    FILE *f;
    f = fopen(name, "w");
    fprintf(f,"%i\n",n);
    for(i=0;i<n;i++) 
        fprintf(f,"%i,%i,%lf\n",I[i],J[i],V[i]);
    fclose(f);
}

int read_csv_size(char* name)
{
    int n;
    FILE *f;
    f = fopen(name, "r");
    fscanf(f,"%i\n",&n);
    fclose(f);
    return n;
}

int read_csv_triplets(char* name, int* I, int* J, double* V)
{
    int n;
    FILE *f;
    f = fopen(name, "r");
    fscanf(f,"%i\n",&n);
    int i = 0;
    for(i=0;i<n;i++)
    {
       fscanf(f,"%i,%i,%lf\n",I+i,J+i,V+i); 
    }
   
   fclose(f);

}

int read_vector(char* name, double* v)
{
    int n;
    FILE *f;
    f = fopen(name, "r");
    fscanf(f,"%i\n",&n);
    int i = 0;
    for(i=0;i<n;i++)
    {
       fscanf(f,"%lf\n",I+i,J+i,V+i); 
    }
   
   fclose(f);
}
