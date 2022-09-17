#include <iostream>
#include <fstream>
using namespace std;

int main() {
    int n = 30;
    int arr[n];
    for (int i = 1; i <= n; i++) {
        arr[i - 1] = 1/i;
    }
    for (int i = 0; i < n; i++) {
        cout << arr[i] << " ";
    }

    return 0;
}