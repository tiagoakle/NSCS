[nf,war] = loadlibrary('libmatlab_shared','matlab_shared.h')
libfunctions('libmatlab_shared') 
%Now try calling get vals
n   = 10;
%Get a pointer to return memory
lpOut = libpointer('doublePtr',zeros(n,1));
%Generate a vector of data
Dat = randn(n,1);
%Call 
res = calllib('libmatlab_shared','CopyVals',lpOut,Dat,n);
%Compare the values
lpOut.Value == Dat

%Test the transpose in triplet function
A = sprand(10,10,0.5);
[I,J,V] = find(A);
lpI = libpointer('int32Ptr',I);
lpJ = libpointer('int32Ptr',J);
nnz = size(I,1);
res = calllib('libmatlab_shared','TransposeInTripletForm',lpI,lpJ,V,nnz);

unloadlibrary libmatlab_shared
