#include <iostream>
#include <cmath>
#include<iomanip>
#include<fstream>
#include<string>

using namespace std;

ofstream fout;
ifstream fin;

string file= "ElectroBoxData";

///constants
    const double pi = 3.14159265;
    const double Ke= 8.98756e9;
    const double Q1 = 1.60217657e-19;
    const double Q2 = 1.60217657e-19;
    const double m = 9.10938291e-31;
    const double k = Ke*Q1*Q2/(m);
    const double k2 = Ke*Q1*Q2;

///to map the particles, initialize a matrix composed of 2 arrays,
/// one for the x coordinate of each particle and one for the y coordinate,
/// to make a square ring that looks like:

///o   o    o    o    o    o    o    o    o
///
///o                                      o
///
///o                                      o
///
///o                                      o
///
///o                                      o
///
///o                                      o
///
///o                                      o
///
///o                                      o
///
///o                                      o
///
///o   o    o    o    o    o    o    o    o

double particleMatrix[32][2]={
                                {10  ,   0},
                                {10  , 2.5},
                                {10  ,   5},
                                {10  , 7.5},
                                {10  ,  10},
                                {7.5 ,  10},
                                {5   ,  10},
                                {2.5 ,  10},
                                {0   ,  10},
                                {-2.5,  10},
                                {-5  ,  10},
                                {-7.5,  10},
                                {-10 ,  10},
                                {-10 , 7.5},
                                {-10 ,   5},
                                {-10 , 2.5},
                                {-10 ,   0},
                                {-10 ,-2.5},
                                {-10 ,  -5},
                                {-10 ,-7.5},
                                {-10 , -10},
                                {-7.5, -10},
                                {  -5, -10},
                                {-2.5, -10},
                                {   0, -10},
                                { 2.5, -10},
                                {   5, -10},
                                { 7.5, -10},
                                {  10, -10},
                                {10  ,-7.5},
                                {10  ,  -5},
                                {10  ,-2.5}
                                            };

double thetaArray [32];
void fThetaArrayGen(double x, double y);

double distancesMatrix[32][2];
void distancesMatrixGen(double x, double y);

double potentialSum();
double potentialGen(double x, double y);
double potentialMagArray[32];
void potentialArrayGen();

double accelSum(bool xy);
double accelGen(double x, double y);
double accelMagArray[32];
void accelMagArrayGen();
///
int main()
{
    fout.open(file.c_str());

    fout <<" i \tx \ty \tv \tKE \tU \tTotalEnergy"<<endl;///print header to file

    double dt=0.000001, x=-0.5, y=0.5, v=25;
    double theta1=-30*pi/180;

    double vx = v*cos(theta1);///calculate x component of initial velocity
    double vy = v*sin(theta1);///calculate y component of initial velocity

    distancesMatrixGen(x,y);///calculate the distances between the electron and each stationary partice
    fThetaArrayGen(x, y);///calculate the angle from each particle to the electron

    double ax=accelSum(1);///sum up the forces caused by each particle on the electron, in the x direction
    double ay=accelSum(0);///sum up the forces caused by each particle on the electron, in the y direction

    double U=potentialSum();///calculate the initial total electric potential energy
    double KE=1/2*m*v*v;///calculate the total initial kinetic energy
    double totalEnergy=U+KE;

    int counter=0;///initialize a counter, for ease in decreasing the amount of data output
    int breaker=0;

for (double t=0;; t+=dt)///loop, adding 1 microsecond to t each iteration
    {
    counter++;

        x += vx*dt+(1/2*ax*dt*dt);///reimman sum for displacement, x direction
        y += vy*dt+(1/2*ay*dt*dt);///reimman sum for displacement, y direction

        vx += ax*dt;///reimman sum for velocity, x direction
        vy += ay*dt;///reimman sum for velocity, y direction
        v=sqrt(vx*vx+vy*vy);

        fThetaArrayGen(x,y);///recalculate the angle from each particle to the electron
        distancesMatrixGen(x,y);///recalculate the distances between the electron and each stationary particle

        ax=accelSum(1);///sum up the forces caused by each particle on the electron at that point, in the x direction
        ay=accelSum(0);///sum up the forces caused by each particle on the electron at that point, in the y direction

        U=potentialSum();///calculate the current total electric potential energy
        KE=(vx*vx+vy*vy)*m*1/2;///calculate the current total kinetic energy
        totalEnergy=U+KE;

        if (counter%1000==0)///every thousandth iteration, to decrease data output
            {
                fout <<t<<"\t"<<x<<"\t"<<y<<"\t"<<v<<"\t"<<KE<<"\t"<<U<<"\t"<< totalEnergy<<endl;

            if (abs(x)>10||abs(y)>10)///if the electron leaves the box, only check every 1000 iterations
            {
                breaker++;///start counting iterations
            }
            }

        if (breaker>40)
            {
                break;/// stop loop once 40 iterations of data have been printed out since the electron left
            }
        }///end for loop
    return 0;
}///end main

void fThetaArrayGen(double x, double y)///recalculated each iteration, updates array containing
                                       ///angle between electron and each particle
{
    for (int i=0;i<32;i++)
    {
        thetaArray[i]=atan(distancesMatrix[0][i]/distancesMatrix[1][i]);
    }
}

double accelGen(double x, double y)///
{
    return k/(x*x+y*y);
}

void distancesMatrixGen(double x, double y)///recalculated each iteration, updates a matrix
                                           ///giving the vector distance between the electron
                                           /// and each particle
{
    for (int i=0; i<32; i++)
    {
        distancesMatrix[i][0]=x-particleMatrix[i][0];
        distancesMatrix[i][1]=y-particleMatrix[i][1];
    }
}
double accelSignGen(double a,double r)///called by accelSum
{

    if (r>0)///if the electron is located to the right of or above the particle
        {
            return a;///a is in the positive direction
        }
    else if (r<0)///if the electron is located to the left of or below the particle
        {
            return -a;///a is in the negative direction
        }
    else///if the electron has the same x or y position as the particle
        {
            return 0;///a is 0
        }
}

double accelSum(bool xy)///recalculated each iteration, finds the net acceleration of electron
{
    accelMagArrayGen();
    double sum=0;
    if (xy==1)
    {
        for (int i=0; i<32; i++)/// go through array and add up all the entries
        {
            sum+=accelSignGen(abs(cos(thetaArray[i])*accelMagArray[i]),distancesMatrix[i][0]);
        }
    }

    else if (xy==0)
    {
        for (int i2=0; i2<32; i2++)/// go through array and add up all the entries
        {
            sum+=accelSignGen(abs(sin(thetaArray[i2])*accelMagArray[i2]),distancesMatrix[i2][1]);
        }
    }
    return sum;
}

void accelMagArrayGen()///called by accelSum, updates an array consisting of the accelerations
                       ///each particle wants to give the electron
{

    for (int i=0; i<32; i++)
    {
        accelMagArray[i]=accelGen(distancesMatrix[i][0],distancesMatrix[i][1]);
    }
}

double potentialSum()///recalculated each iteration, finds total potential energy of electron
{
    potentialArrayGen();
    double sum=0;
        for (int i=0; i<32; i++)/// go through array and add up all the entries
        {
            sum+=potentialMagArray[i];
        }
    return sum;
}

void potentialArrayGen()///called by potentialSum, updates an array consisting of U
                        ///of the electron due to each point particle
{

    for (int i=0; i<32; i++)
    {
        potentialMagArray[i]=potentialGen(distancesMatrix[i][0],distancesMatrix[i][1]);
    }
}

double potentialGen(double x, double y)///called by potentialArrayGen, returns magnitude of U
{
    return k2/sqrt(x*x+y*y);
}
