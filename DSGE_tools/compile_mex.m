%*******************************************************************
% COMPILE MEX FILES - WORKS ON A MAC OSX 10.9 FOR SURE
%
% Uncomment/comment to compile only particular files
%*******************************************************************
mex -v -largeArrayDims dim5_simplex_eval_mex.c  -lmwlapack -lmwblas

mex -v -largeArrayDims dim6_simplex_eval_mex.c  -lmwlapack -lmwblas

mex -v -largeArrayDims dim7_simplex_eval_mex.c  -lmwlapack -lmwblas

mex -v -largeArrayDims dim8_simplex_eval_mex.c  -lmwlapack -lmwblas



%mex -v -largeArrayDims uf_c_mex.c  -lmwlapack -lmwblas
%mex -v -largeArrayDims test_sparse_mex.c  -lmwlapack -lmwblas