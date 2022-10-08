#include <iostream>
#include <fstream>
using namespace std;

// parameters of simulation
int n = 10000;
double dt = 0.01;
//--------------------------

class Object {
public:
    double x_0 = 10, v_0 = 0, w = 1;
    double* x = new double[n];
    double* v = new double[n];
    double* energy = new double[n];
    Object(){
        *x = x_0;
        *v = v_0;}
    void measureEnergy() {
        for (int i = 0; i < n; i++) {
            energy[i] = v[i] * v[i] / 2 + w * w * x[i] * x[i] / 2;
        }
    }

};

class ModMethod{
public:
    Object object;

public:
    void EulersMethod(){
        for (int i = 0; i < n; i++){
            object.x[i + 1] = object.x[i] + dt * object.v [i];
            object.v[i + 1] = object.v[i] + dt * (-object.w*object.w) * object.x[i];
        }
    }

    void HoynsScheme(){
        auto* x_1 = new double[n];
        auto* v_1 = new double[n];
        for (int i = 0; i < n; i++){
            x_1[i + 1] = x_1[i] + dt * v_1[i];
            v_1[i + 1] = v_1[i] + dt * (-object.w * object.w) * x_1[i];
        }
        for (int i = 0; i < n; i++) {
            object.x[i + 1] = object.x[i] + dt/2 * (object.v[i] + v_1[i + 1]);
            object.v[i + 1] = object.v[i] + dt/2 * (-object.w * object.w)*(object.x[i] + x_1[i + 1]);
        }
    }

};

class FileOutput{
public:
    Object object;
public:
    void write(){
        ofstream out;
        out.open("/home/starman/CLionProjects/EulersMethod/data.txt");
        if (out.is_open())
        {
            for (int i = 0; i < n; i++){
                out <<  dt * i  << ' ' << object.x[i]<< ' ' << object.energy[i] << ' ' << object.v[i] <<'\n';
            }
        }
        out.close();
    }
};

int main() {
    Object object;
    ModMethod method;
    method.object = object;
    method.HoynsScheme();
    object.measureEnergy();
    FileOutput file;
    file.object = object;
    file.write();

    return 0;
}
