#include <matrix.h>
#include <mex.h>

void MyOp_conj_trans_z_rot(double *blur, double *orig, double *index,double *weight, int length)
{
	int i,j,ind;
	for(i=0;i<length;i++)
	{
        j=index[i]-1;
        if(j>=0 && j<length)
        {
            blur[j]=blur[j]+orig[i]*weight[i];
        }
        ind=i+length;
		j=index[ind]-1;
        if(j>=0 && j<length)
        {
            blur[j]=blur[j]+orig[i]*weight[ind];
        }
		ind=i+2*length;
		j=index[ind]-1;
        if(j>=0 && j<length)
        {
            blur[j]=blur[j]+orig[i]*weight[ind];
        }
		ind=i+3*length;
		j=index[ind]-1;
        if(j>=0 && j<length)
        {
            blur[j]=blur[j]+orig[i]*weight[ind];
        }
	}
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	double *blur,*orig,*weight,*index;
	int mrows0,mcols0,mrows1,mcols1,mrows2,mcols2;
	if(nrhs<2)
	{
		mexErrMsgTxt("Require at least three inputs");
	}
	mrows0=mxGetM(prhs[0]);
	mcols0=mxGetN(prhs[0]);
	mrows1=mxGetM(prhs[1]);
	mcols1=mxGetN(prhs[1]);
	mrows2=mxGetM(prhs[2]);
	mcols2=mxGetN(prhs[2]);
	if(mcols0!=mcols1/4 || mcols1!=mcols2 || mrows0!=1 || mrows1!=1 || mrows2!=1)
		mexErrMsgTxt("Please input three 1*n 1*4n 1*4n vectors");

	plhs[0]=mxCreateDoubleMatrix(mrows0,mcols0,mxREAL);

	blur=mxGetPr(plhs[0]);
	orig=mxGetPr(prhs[0]);
	index=mxGetPr(prhs[1]);
	weight=mxGetPr(prhs[2]);

	MyOp_conj_trans_z_rot(blur,orig,index,weight,mcols0);
}