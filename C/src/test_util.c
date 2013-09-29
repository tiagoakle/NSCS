#include <stdio.h>
#include "test_util.h"

int read_csv_size(char* name,int* m, int *n, int* nnz)
{
    
    FILE *f;
    f = fopen(name, "r");
    fscanf(f,"%i,%i,%i\n",m,n,nnz);
    fclose(f);
    return 0;
}

int read_csv_triplets(char* name, int* I, int* J, double* V)
{
    int n;
    int m;
    int nnz;

    FILE *f;
    f = fopen(name, "r");
    fscanf(f,"%i,%i,%i\n",&m,&n,&nnz);
    
    int i = 0;
    for(i=0;i<nnz;i++)
    {
       fscanf(f,"%i,%i,%lf\n",I+i,J+i,V+i); 
    }
   
   fclose(f);

}

int read_size(char* name)
{
    int n;
    FILE *f;
    f = fopen(name, "r");
    fscanf(f,"%i\n",&n);
    fclose(f);
    return n;
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
       fscanf(f,"%lf\n",v+i); 
    }
   
   fclose(f);
}

int write_matrix_to_csv(char* name, int* I, int* J, double* V, int nnz, int m, int n)
{
    int i=0;
    FILE *f;
    f = fopen(name, "w");
    fprintf(f,"%i,%i,%i\n",m,n,nnz);
    for(i=0;i<nnz;i++) 
        fprintf(f,"%i,%i,%lf\n",I[i],J[i],V[i]);
    fclose(f);
}

int write_vector_to_csv(char* name, double* v, int n)
{
    int i=0;
    FILE *f;
    f = fopen(name, "w");
    fprintf(f,"%i\n",n);
    for(i=0;i<n;i++) 
    {
        fprintf(f,"%lf\n",v[i]);
    }
    fclose(f);
}

//Reads the vectors written by the matlab write_vector_bin.m file
int read_vector_bin_size(char* fname)
{
    int size; 
    FILE* f = fopen(fname,"rb");
    if(f)
    {
        fseek(f,0,SEEK_END);
        size = ftell(f);
    }
    fclose(f);
    return size/sizeof(double); 
     
}

void read_vector_bin(char* fname, double* x, int n)
{
    int size; 
    FILE* f = fopen(fname,"rb");
    if(f)
    {
        fread((void*)x,sizeof(double),n,f);
    }else 
        printf("Unable to open %s \n",fname); 
    fclose(f);
}

//Reads the matrices written by the matlab write_vector_bin.m file
void read_matrix_bin_size(char* fname, int* m, int* n, int* nnz)
{
    int size; 
    FILE* f = fopen(fname,"rb");
    int read; 
    if(f)
    {
       read =  fread((void*)m,sizeof(int),1,f);
       if(read!=1){fprintf(stderr,"Unable to read from %s\n",fname);}
       read =  fread((void*)n,sizeof(int),1,f);
       if(read!=1){fprintf(stderr,"Unable to read from %s\n",fname);}
       read =  fread((void*)nnz,sizeof(int),1,f);
       if(read!=1){fprintf(stderr,"Unable to read from %s\n",fname);}
    }else 
        printf("Unable to open %s \n",fname); 
    fclose(f);
     
}

void read_matrix_bin(char* fname, int* Ai, int* Aj, double* Av)
{
    int m,n,nnz;
    FILE* f = fopen(fname,"rb");
    if(f)
    {
        fread((void*)&m,sizeof(int),1,f);
        fread((void*)&n,sizeof(int),1,f);
        fread((void*)&nnz,sizeof(int),1,f);
        fread((void*)Ai,sizeof(int),nnz,f);
        fread((void*)Aj,sizeof(int),nnz,f);
        fread((void*)Av,sizeof(double),nnz,f);
    } else 
        printf("Unable to open %s \n",fname);

    fclose(f); 
}

