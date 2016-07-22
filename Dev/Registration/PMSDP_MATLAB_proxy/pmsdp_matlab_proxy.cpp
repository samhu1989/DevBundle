#include "pmsdp_matlab_proxy.h"
#include "libPMSDP.h"
bool initMatlab(void)
{
    if(!libPMSDPInitialize())
    {
        std::cerr<<"Failed to Initialize Matlab C Runtime"<<std::endl;
        return false;
    }
    return true;
}
void compute(
        const arma::mat &P,
        const arma::mat &Q,
        arma::mat& R,
        arma::uvec& X
        )
{
    assert( P.n_rows == Q.n_rows );
    mxArray* plhs[2];
    mxArray* prhs[3];
    prhs[0] = mxCreateDoubleMatrix(P.n_rows,P.n_cols,mxREAL);
    memcpy(mxGetPr(prhs[0]),(mxDouble*)P.memptr(),sizeof(double)*P.size());
    prhs[1] = mxCreateDoubleMatrix(Q.n_rows,Q.n_cols,mxREAL);
    memcpy(mxGetPr(prhs[1]),(mxDouble*)Q.memptr(),sizeof(double)*Q.size());
    /*
    params.probDim = probDim;
    params.n = n;
    params.k = k;
    params.UtilizeXflag = UtilizeXflag;
    params.permConstraint = permConstraint;
    params.utilizeRFlag = utilizeRFlag;
    params.Rtol = Rtol;
    params.verbose = verbose;
    */
    size_t ndim = 1, dims[1] = {1};
    size_t number_of_fields = 8;
    const char *field_names[] = {
        "probDim",
        "n",
        "k",
        "UtilizeXflag",
        "permConstraint",
        "utilizeRFlag",
        "Rtol",
        "verbose"
    };
    prhs[2]=mxCreateStructArray(ndim,dims,number_of_fields,field_names);
    mxSetField(prhs[2],0,field_names[0],mxCreateDoubleScalar(P.n_rows));
    mxSetField(prhs[2],0,field_names[1],mxCreateDoubleScalar(P.n_cols));
    mxSetField(prhs[2],0,field_names[2],mxCreateDoubleScalar(Q.n_cols));
    mxSetField(prhs[2],0,field_names[3],mxCreateLogicalScalar(false));
    mxSetField(prhs[2],0,field_names[4],mxCreateDoubleMatrix(0,0,mxREAL));
    mxSetField(prhs[2],0,field_names[5],mxCreateLogicalScalar(false));
    mxSetField(prhs[2],0,field_names[6],mxCreateDoubleScalar(-1.0));
    mxSetField(prhs[2],0,field_names[7],mxCreateLogicalScalar(false));
    mlxSolvePMSDP(2,plhs,3,prhs);
    R = arma::mat((double*)mxGetData(plhs[1]),mxGetM(plhs[1]),mxGetN(plhs[1]),true,true);
    arma::mat dX((double*)mxGetData(plhs[0]),mxGetM(plhs[0]),mxGetN(plhs[0]),false,true);
    X = arma::uvec(dX.n_cols);
    #pragma omp parallel for
    for(arma::uword ic=0;ic < dX.n_cols;++ic)
    {
        arma::uvec index = arma::find(dX.col(ic)==0);
        if(!index.is_empty())X(ic) = index(0);
    }
    mxDestroyArray(prhs[0]);
    mxDestroyArray(prhs[1]);
    mxDestroyArray(prhs[2]);
    mxDestroyArray(plhs[0]);
    mxDestroyArray(plhs[1]);
}
void terminateMatlab(void)
{
    libPMSDPTerminate();
}
