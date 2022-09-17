#include <iostream>
#include <fstream>
using namespace std;

int fibbonacci(int num){
    if (num == 1 || num == 2) return 1;
    return fibbonacci(num-1)+fibbonacci(num-2);
}

int main() {
    int n = 30;
    int arr[n];
    for (int i = 1; i <= n; i++){
        int fib = fibbonacci(i);
        arr[i-1] = fib;
    }
    for (int i = 0; i < n; i++){
        cout << arr[i] << " ";
    }

    return 0;
}