#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<stdint.h>

int16_t multFixed(int16_t val1, short int q1, int16_t val2, short int q2, short int format){
    int16_t res;
    res = ((int) val1 * val2)/(1<<(q1 + q2 - format));
    return res;
}

int main(){
    float val1 = 0.03125; short int q1 = 12;
    float val2 = 0.03125; short int q2 = 11;
    short int finalFormat = 11;
    float resultFloat;

    int16_t fixedVal1 = (int16_t)(val1 * (1<<q1));
    int16_t fixedVal2 = (int16_t)(val2 * (1<<q2));
    int16_t result;

    result = multFixed(fixedVal1, q1, fixedVal2, q2, finalFormat);

    resultFloat = (float)result / (1<<finalFormat);
    printf("%f * %f = %f", val1, val2, resultFloat);

    return 0;
}