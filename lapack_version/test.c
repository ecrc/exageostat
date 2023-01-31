#include <vector>
#include <algorithm>
#include<iostream>
#include<cmath>
#include "time.h"
using namespace std;

int main () 
{
   
clock_t begin = clock();
 int n=400000,  m=1000;  
    double x=0,y=0;
    double s=0;
    vector< double > shifts(n,0);


    #pragma omp parallel for 
    for (int j=0; j<n; j++) {

        double r=0.0;
        for (int i=0; i < m; i++){

            double rand_g1 = m;
            double rand_g2 = m;     

            x += rand_g1;
            y += rand_g2;
            r += sqrt(rand_g1*rand_g1 + rand_g2*rand_g2);
        }
        shifts[j] = r / m;
    }

    cout << *std::max_element( shifts.begin(), shifts.end() ) << endl;



/* here, do your time-consuming job */

clock_t end = clock();
double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
cout<<"time spent: "<<time_spent<<endl;

}
