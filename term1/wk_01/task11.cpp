#include <iostream>
#include <fstream>
using namespace std;

int main() {
    ofstream file("D:\\task11.bin", ios::out | ios::binary);
   if (file.is_open())
   {
       int x = 0x42ABCDEF;
       cout << "x = " << x << "\n";
       file.write((char *) &x, sizeof(x));
       file.close();
   };

   return 0;
}