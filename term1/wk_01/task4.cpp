#include <iostream>
#include <fstream>
using namespace std;
int main() {
    ofstream out;
    out.open("D:\\task4.txt");
    if (out.is_open())
    {
        for (int i = 1; i < 31; i++){
            out << i << " ";
        }
    }
    out.close();
    return 0;
}

