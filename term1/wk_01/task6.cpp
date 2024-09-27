#include <iostream>
#include <fstream>
#include <string>
using namespace std;
int main(int argc, char *argv[]) {

    if (argc != 2){
        std::cerr << "Error: wrong number of arguments" << std::endl;
        return -1;
    }

    int n = stoi(argv[1]);

    ofstream out;
    out.open("D:\\task6.txt");
    if (out.is_open())
    {
        for (int i = 1; i <= n; i++){
            out << i << " ";
        }
    }
    out.close();
    return 0;
}