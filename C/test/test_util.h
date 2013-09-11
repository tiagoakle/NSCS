int read_csv_size(char* name,int* m, int *n, int* nnz);
int read_csv_triplets(char* name, int* I, int* J, double* V);
int read_vector(char* name, double* v);
int read_size(char* name);
int write_matrix_to_csv(char* name, int* I, int* J, double* V, int nnz, int m, int n);
int write_vector_to_csv(char* name, double* v, int n);
int read_vector_bin_size(char* fname);
void read_vector_bin(char* fname, double* x, int n);
