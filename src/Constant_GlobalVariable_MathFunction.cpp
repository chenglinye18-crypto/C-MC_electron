#include <stdio.h>
#include "Constant_GlobalVariable_MathFunction.h"
#include <sstream>

double sia0,sirho,siul,siut,a0pi,sieg;

double T0,eV0,em0,hq0,ec0,rmom0,spr0,spk0,time0,velo0, Nc, Ncc, psi_si, phi_top, N_cur;

double mx[3];
double my[3];
double mz[3];
double Nccc[3];
double ml, mt;

double cvr,pot0,field0,conc0,dens0,dpc0,scrt0,curr0;

double ChargeSign[NumParType]={-1.0,1.0,-1.0};

bool gaasfl;

bool sifl;
bool gefl;


double MASSeff[NumParType]={0.286,0.8,1.0};

double Max(double p1,double p2)
{
    double a;
    a=p1>p2 ? p1:p2;
    return a;
}


 double Max(double p1,double p2,double p3)
{
    double a;
    a=p1>p2?p1:p2;
    a=a>p3?a:p3;
    return a;
}


 double Max(double p1,double p2,double p3,double p4)
{
    double a;
    a=p1>p2?p1:p2;
    a=a>p3?a:p3;
    a=a>p4?a:p4;
    return a;
}


 double Max(double p1,double p2,double p3,double p4,double p5)
{
    double a;
    a=p1>p2?p1:p2;
    a=a>p3?a:p3;
    a=a>p4?a:p4;
    a=a>p5?a:p5;
    return a;
}


 double Min(double p1,double p2)
{
    double a;
    a=p1<p2?p1:p2;
    return a;
}


 double Min(double p1,double p2,double p3)
{
    double a;
    a=p1<p2?p1:p2;
    a=a<p3?a:p3;
    return a;
}


 double Min(double p1,double p2,double p3,double p4)
{
    double a;
    a=p1<p2?p1:p2;
    a=a<p3?a:p3;
    a=a<p4?a:p4;
    return a;
}


 double Min(double p1,double p2,double p3,double p4,double p5,double p6, int & flag)
{
    double a;
    a = p1; 
    flag = 1;
    if (a > p2) {
      a = p2;
      flag = 2;
    }
    if (a > p3) {
      a = p3;
      flag = 3;
    }
    if (a > p4) {
      a = p4;
      flag = 4;
    }
    if (a > p5) {
      a = p5;
      flag = 5;
    }
    if (a > p6) {
      a = p6;
      flag = 6;
    }
    return a;
}

