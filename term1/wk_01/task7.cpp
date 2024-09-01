#include <iostream>
#include <fstream>
using namespace std;

int fibbonacci(int num){
    if (num == 1 || num == 2) return 1;
    return fibbonacci(num-1)+fibbonacci(num-2);
}

int main() {
    int n = 30;
    ofstream out;
    out.open("D:\\task7.txt");
    if (out.is_open())
    {
        for (int i = 1; i <= n; i++){
            out << i << "      " << fibbonacci(i) << endl;
        }
    }
    out.close();
    return 0;
}