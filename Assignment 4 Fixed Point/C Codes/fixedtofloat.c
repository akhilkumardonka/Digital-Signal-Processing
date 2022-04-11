#include<stdlib.h>
#include<stdio.h>
#include<stdint.h>
#include<string.h>

int main(){
    int16_t fxpValue = -462;
    float floatValue;
    int16_t qValue = 7;

    floatValue = (float)fxpValue/(1<<qValue);
    printf("float value for fixed point = %d with q = %d is : %f", fxpValue, qValue, floatValue);

    return 0;
}