#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<stdint.h>

int main(){
    int16_t fxpValue;
    float floatValue = -3.613;
    uint8_t qValue = 7;

    fxpValue = (int16_t)(floatValue*(1<<qValue));
    printf("fixed point value for %f, with Q = %d is : %d",floatValue, qValue, fxpValue);

    return 0;
}

