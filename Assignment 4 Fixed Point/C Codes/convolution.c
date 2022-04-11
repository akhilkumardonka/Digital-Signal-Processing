#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define MAXLEN 100

int fxdAddScalar(int num1, int num2)
{
    return num1+num2;
}

int bitshift(int num,int shift)
{
    return num*(1<<shift);
}

int fxdMulScalar(int num1, int q1, int num2, int q2, int resq)
{
    long long temp = (long long) num1*num2;
    int res;
    res = temp/(1<<(q1+q2-resq));
    return res;
}

float *convolutionflt(float X[MAXLEN], float H[MAXLEN], float Y[MAXLEN], int lenX, int lenH, int lenY){
    for(int n=0;n<lenY;n++)
    {
        Y[n] = 0;
        for(int k = 0;k<=n;k++)
        {
        if((n-k)>=lenH || k>=lenX)
        continue;
        Y[n] = Y[n] + (X[k]*H[n-k]);
        }
    }
    return Y;
}

int *convolutionfxd(int X[MAXLEN], int H[MAXLEN], int Y[MAXLEN], int lenX, int lenH, int lenY, int Q)
{
    for(int n=0;n<lenY;n++)
    {
        Y[n] = 0;
        for(int k = 0;k<=n;k++)
        {
        if((n-k)>=lenH || k>=lenX)
        continue;
        Y[n] = fxdAddScalar(bitshift(Y[n],0),fxdMulScalar(X[k],Q,H[n-k],Q,Q));
        }
    }
    return Y;
}

int main(){
    float X[MAXLEN];
    float H[MAXLEN];
    float Y[MAXLEN];

    int Xfix[MAXLEN];
    int Hfix[MAXLEN];
    int Yfix[MAXLEN];

    int xlen = 10;
    int hlen = 25;
    int ylen = xlen +hlen -1;
    float *floatyConv;
    int *fixedyConv;
    float fixedToFloatY[MAXLEN];

    int niter = 10000;
    int qValue[] = {7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27.28,29,30,31};
    float errors[MAXLEN];

    for(int q = 0;q<sizeof(qValue)/sizeof(qValue[0]);q++){

        float iterationError=0;

        for(int j =0; j<niter;j++){

            float err=0;

            // Generating random numbers for X
            for (int i=0;i<xlen;i++)
            {
                X[i] = ((float)rand()/RAND_MAX)*(float)(10);
            }

            // Generating random numbers for H
            for(int i=0;i<hlen;i++)
            {
                H[i] =  ((float)rand()/RAND_MAX)*(float)(10);
            }

            // Fixed point conversion for X
            for (int i=0;i<xlen;i++)
            {
                Xfix[i] = (int)(X[i]*(1<<qValue[q]));
            }

            // Fixed point conversion for H
            for (int i=0;i<hlen;i++)
            {
                Hfix[i] = (int)(H[i]*(1<<qValue[q]));
            }
            
            // floating point convolution
            floatyConv=convolutionflt(X, H, Y, xlen, hlen, ylen);

            // fixed point convolution
            fixedyConv=convolutionfxd(Xfix, Hfix, Yfix, xlen, hlen, ylen, qValue[q]);

            // fixed point convolution output to floating point values
            for(int i=0; i<ylen; i++){
                fixedToFloatY[i] = (float)fixedyConv[i]/(1<<qValue[q]);
                // printf("Float Y: %f | Fixed Y: %d | Fixed to Float Y: %f\n", floatyConv[i], fixedyConv[i], fixedToFloatY[i]);
            }

            for(int i=0; i<ylen; i++){
                err += pow((floatyConv[i] - fixedToFloatY[i]), 2);
            }

            err = err/ylen;
            iterationError += err;
        }
        
        errors[q] = iterationError/niter;
        printf("Q : %d, Error : %f\n", qValue[q], iterationError/niter);
    } 

    

    return 0;
}