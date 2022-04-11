#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<stdint.h>

int main(){
    float val1 = 0.03125; short int q1 = 12;
    float val2 = 0.03125; short int q2 = 11;
    short int format;

    int16_t fxpVal1 = (int16_t)(val1*(1<<q1));
    int16_t fxpVal2 = (int16_t)(val2*(1<<q2));

    if(q1 > q2){
        fxpVal2 = fxpVal2*(1<<(q1-q2));
        format = q1;
    }else{
        fxpVal1 = fxpVal1*(1<<(q2-q1));
        format = q2;
    }

    int16_t result = fxpVal1 + fxpVal2;

    float resultFloat = (float)result / (1<<format);
    printf("Value 1 : %f \n", val1);
    printf("Value 2 : %f \n", val2);
    printf("Value 1 + Value 2 = %f", resultFloat);

    return 0;
}