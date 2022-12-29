#include <iostream>
#include "mat.h"

using namespace std;

int main(){
    matriz <double> A(3,3,{
        0,2,3,
        4,0,0,
        7,8,9
    });
    cout << A.det() << endl;
    return 0;
}