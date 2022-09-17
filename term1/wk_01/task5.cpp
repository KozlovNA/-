#include <iostream>
#include <fstream>
using namespace std;
int main() {
    int n;
    cin >> n;
    ofstream out;
    out.open("D:\\task5.txt");
    if (out.is_open())
    {
        for (int i = 1; i <= n; i++){
            out << i << " ";
        }
    }
    out.close();
    return 0;
}