#include <iostream>
#include "mat.h"

using namespace std;

int main(){
    matriz <double> A(3,3,{
        0,2,3,
        4,0,0,
        7,8,9
    });

    std::cout << std::setprecision(5);
    std::cout << A.inv() << endl;
    cout << A.norma() << endl;
    return 0;
}