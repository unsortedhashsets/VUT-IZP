/*******************************************************
* Soubor:  proj2.c                                     *
* Datum:   19.11.2018                                  *
* Autor:   Abramov Mikhail, xabram00@stud.fit.vutbr.cz *
* Projekt: IZP c. 2 - Iteracni vypocty                 *
*******************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
void printHelp(void);
double taylor_log(double x, unsigned int n);
double cfrac_log(double x, unsigned int n);
double taylor_pow(double x, double y, unsigned int n);
double taylorcf_pow(double x, double y, unsigned int n);
void printHelp(void) {
    printf("Iterative calculations program\n"
           "By Abramov Mikhail\n"
           "The program performs the calculation of the entered mathematical function\n"
           "with certain number of iterations. Program reads characters from standard input.\n"
           "Prints results on the screen.\n"
           "Use: proj2 --h\n"
           "     proj2 --log x n|-pow x y n\n"
           "Popis parametru:\n"
           "  --h           Display this help screen.\n"
           "  --log x n     Calculates the natural logarithm for x by the number of iterations n\n"
           "  --pow x y n   Calculates x to the power of Y with the number of iterations n\n");
}
double taylorcf_pow(double x, double y, unsigned int n){
    if (x==0){//X validation (if x=0 result = 0)
        if (y==0) { //X validation (if x=0 y =0 result = 1)
            return 1;
        }
        if (y<0){
            return INFINITY; // because of 1/0
        }
        return 0;
    }
    if (x<0 && fmod(y,1.0)!=0){//X validation (if x<0 a Y is decimal fractions (-3^2.5)
        return -NAN;
    }
    int SX=0;//FLAG for +-x
    int SY=0;//FLAG for +-y
    if (x<0){// if x<0 change the flag to 1 and y to +
        x=-x;
        SX=1;
    }
    if (y<0){// if y<0 change the flag to 1 and x to +
        y=-y;
        SY=1;
    }
    double lnX=(cfrac_log(x,n));//calculation of ln in cf
    double yNew=1;//Realisation of Tlp for POW
    double yOld=1;
    for (unsigned int i=1; i<n; i++) {
        yOld*=y*lnX/i;
        yNew+=yOld;
    }
    if (SX!=0){//validation of FLAGS and +-X +-Y
        if (SY!=0 && fmod(y,2.0)==0.0 ) {//if x-yFlAGS =1 and y is whole => result will be 1/result
            yNew = 1/yNew;
            return yNew;
        }
        else if (SY!=0){//if x-yFlAGS =1 result will be -1/result
            yNew=1/yNew;
            yNew=-yNew;
            return yNew;
        }
        else if (fmod(y,2.0)==0.0 ){//if x-FlAG =1 and y is whole result will be result
            yNew=yNew;
            return yNew;
        }
        yNew=-yNew;
        return yNew;
    }
    else {
        if (SY!=0) {
            yNew=1/yNew;
            return yNew;
        }
    return yNew;
    }
}
double taylor_pow(double x, double y, unsigned int n) {
    if (x==0){//X validation (if x=0 result = 0)
        if (y==0) { //X validation (if x=0 y =0 result = 1)
            return 1;
        }
        if (y<0){
            return INFINITY; // because of 1/0
        }
        return 0;
    }
    if (x<0 && fmod(y,1.0)!=0){//X validation (if x<0 a Y is decimal fractions (-3^2.5)
        return -NAN;
    }
    int SX=0;//FLAG for +-x
    int SY=0;//FLAG for +-y
    if (y<0){// if y<0 change the flag to 1 and y to +
        y=-y;
        SY=1;
    }
    if (x<0){// if x<0 change the flag to 1 and x to +
        x=-x;
        SX=1;
    }
    double lnX = (taylor_log(x,n));//calculation of ln in tl
    double yNew=1;//Realisation of Tlp for POW
    double yOld=yNew;
    for (unsigned int i=1; i<n; i++)
    {
        yOld*=y*lnX/i;
        yNew+=yOld;
    }
    if (SX!=0){//validation of FLAGS and +-X +-Y
        if (SY!=0 && fmod(y,2.0)==0.0 ) {//if x-yFlAGS =1 and y is whole => result will be 1/result
            yNew = 1/yNew;
            return yNew;
        }
        else if (SY!=0){//if x-yFlAGS =1 result will be -1/result
            yNew=1/yNew;
            yNew=-yNew;
            return yNew;
        }
        else if (fmod(y,2.0)==0.0 ){//if x-FlAG =1 and y is whole result will be result
            yNew=yNew;
            return yNew;
        }
        yNew=-yNew;
        return yNew;
    }
    else {
        if (SY!=0) {
            yNew=1/yNew;
            return yNew;
        }
        return yNew;
    }
}
double taylor_log(double x, unsigned int n){
    double yNew, yOld;
    if (x==0){
        return -INFINITY;//X validation (if X=0 => result -inf)
    }
    else if (x<1){//REALISATION of TL for x<0
        double y=1-x;
        yNew=-y;
        yOld=y;
        for (unsigned int i=2; i<=n; i++) {
            yOld*=y/i*(i-1);
            yNew-=yOld;
        }
        return yNew;
    }
    else {//REALISATION of TL for x>=0
        yNew=(x-1)/x;
        yOld=yNew;
        for (unsigned int i=2; i<=n; i++) {
            yOld*=((x-1)/x)/i*(i-1);
            yNew+=yOld;
        }
        return yNew;
    }
}
double cfrac_log(double x, unsigned int n){
    if (x==0){//X validation (if X=0 => result -inf)
        return -INFINITY;
    }
    double z=(x-1)/(1+x);//REALISATION OF CF
    double yOld=2*n-1;
    double yNew=2*z;
    for (unsigned int i=n; i>1; i--){
        yOld=(2*i-3)-(z*z*(i-1)*(i-1))/yOld;
    }
    yNew/=yOld;
    return yNew;
}
int main(int argc, char * argv[]) {
    double x, y;
    unsigned int n;
    char *err;
    if (strcmp(argv[1], "--log") == 0) {//first argument validation (function LOG validation)
        if (argc != 4) {
            printf("ERROR (1): the number of function arguments for logarithm is 2:\n"
                   "(y-decimal fraction, N-whole number)\nTo enter instructions --h\n");
            return 0;
        }
        x = strtod(argv[2], &err);//Second argument validation (X validation must be >0 for log)
        if (*err || x < 0) {
            printf("ERROR (2): x - must be decimal fraction greater than or equal ZERO\nTo enter instructions --h\n");
            return 0;
        }
        n = strtoul(argv[3], &err, 10);//Third argument scanning and validation (N validation must be >0 and whole)
        if (*err || n <= 0) {
            printf("ERROR (3): n - must be whole number\nTo enter instructions --h\n");
            return 0;
        }
        printf("       log(X) = %.12g\n", log(x));// X3 results for LOG
        printf(" cfrac_log(X) = %.12g\n", cfrac_log(x, n));
        printf("taylor_log(X) = %.12g\n", taylor_log(x, n));
    } else if (strcmp(argv[1], "--pow") == 0) {//first argument validation (function POW validation)
        if (argc != 5) {
            printf("ERROR (4): the number of function arguments for POW is 3:\n"
                   "(x-decimal fraction, Y-decimal fraction, N-whole number)\nTo enter instructions --h\n");
            return 0;
        }
        x = strtod(argv[2], &err);//Second argument scanning and validation
        if (*err) {
            printf("ERROR (5): x - must be decimal fraction\nTo enter instructions --h\n");
            return 0;
        }
        y = strtod(argv[3], &err);//Third argument scanning and validation
        if (*err) {
            printf("ERROR (5): x - must be decimal fraction\nTo enter instructions --h\n");
            return 0;
        }
        n = strtoul(argv[4], &err, 10);//Forth argument scanning and validation (N validation must be >0 and whole)
        if (*err) {
            printf("ERROR (7): n - must be whole number\nTo enter instructions --h\n");
            return 0;
        }
        printf("          pow(X,Y)= %.12g\n", pow(x, y)); // X3 results for POW
        printf("  taylor_pow(X,Y) = %.12g\n", taylor_pow(x, y, n));
        printf("taylorcf_pow(X,Y) = %.12g\n", taylorcf_pow(x, y, n));
        return 0;
    } else if (strcmp(argv[1], "--h") == 0) {
        {//first argument validation (function HELP validation)
            printHelp();
            return 0;
        }
    } else {//first argument validation (ERROR)
        printf("ERROR (8): Bad Request.\nTo enter instructions --h\n");
        return 0;
    }
}
