#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

// parameters of simulation
double duration = 10;
double dt = 0.01;
int n = floor(duration/dt);
//--------------------------

class Object {
public:
    double x_0, v_0, w;
    double* x;
    double* v;
    double* energy;
    Object():
        x(new double[n]),
        v(new double[n]),
        energy(new double[n]),
        x_0(10), v_0(0), w(1)
        {
            *x = x_0;
            *v = v_0;}

    ~Object()
    {
        delete[] x;
        delete[] v;
        delete[] energy;
    }

    Object& operator=(Object const &src) {

        double *new_x = new double[n];
        for (size_t pos = 0; pos != n; ++pos)
            new_x[pos] = x[pos];
        delete[] x;
        x = new_x;

        double *new_v = new double[n];
        for (size_t pos = 0; pos != n; ++pos)
            new_v[pos] = v[pos];
        delete[] v;
        v = new_v;

        double *new_energy = new double[n];
        for (size_t pos = 0; pos != n; ++pos)
            new_energy[pos] = energy[pos];
        delete[] energy;
        energy = new_energy;

        return *this;
    }

    Object(Object const &src): Object(){
        for (size_t pos = 0; pos != n; pos++)
            x[pos] = src.x[pos];

        for (size_t pos = 0; pos != n; pos++)
            v[pos] = src.v[pos];

        for (size_t pos = 0; pos != n; pos++)
            energy[pos] = src.energy[pos];
    }

    void measureEnergy() {
        for (int i = 0; i < n; i++) {
            energy[i] = v[i] * v[i] / 2 + w * w * x[i] * x[i] / 2;
        }
    }

};

double KahanSum(double* input, int inputsize) {
    double sum = 0.0;
    double c = 0.0;
    for(int i = 1; i < inputsize; i++){
        double y = input[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

double KahanSum(double input, double term, double c) {
        double y = term - c;
        double t = input + y;
        c = (t - input) - y;
    return t;
}


class HoynsSchemeKahan{
public:
    double* x_1;
    double* v_1;

public:
    HoynsSchemeKahan():
        x_1(new double[n]),
        v_1(new double[n]) {}

    ~HoynsSchemeKahan() {delete[] x_1; delete[] v_1;}

void count(Object &object){
    double c_x = 0.0;
    double c_v = 0.0;
    for (int i = 0; i < n - 1; i++){
        x_1[i + 1] = x_1[i] + dt * v_1[i] - c_x;
        c_x = (x_1[i + 1] - x_1[i]) - (dt * v_1[i] - c_x);
        v_1[i + 1] = v_1[i] + dt * (-object.w * object.w) * x_1[i] - c_v;
        c_v = (v_1[i + 1] - v_1[i]) - (dt * (-object.w * object.w) * x_1[i] - c_v);
    }
    c_x = 0.0;
    c_v = 0.0;
    for (int i = 0; i < n - 1; i++) {
        object.x[i + 1] = object.x[i] + dt/2 * (object.v[i] + v_1[i + 1]) - c_x;
        c_x = (object.x[i + 1] - object.x[i]) - (dt/2 * (object.v[i] + v_1[i + 1]) - c_x);
        object.v[i + 1] = object.v[i] + dt/2 * (-object.w * object.w)*(object.x[i] + x_1[i + 1]) - c_v;
        c_v = (object.v[i + 1] - object.v[i]) - (dt/2 * (-object.w * object.w)*(object.x[i] + x_1[i + 1]) - c_v);
    }
}
     HoynsSchemeKahan& operator=(HoynsSchemeKahan const &src) = delete;
     HoynsSchemeKahan(HoynsSchemeKahan const &src) = delete;
};

class HoynsScheme{
public:
    double* x_1;
    double* v_1;

public:
    HoynsScheme():
            x_1(new double[n]),
            v_1(new double[n]) { }

    ~HoynsScheme() {delete[] x_1; delete[] v_1;}

    void count(Object &object){
        for (int i = 0; i < n - 1; i++){
            x_1[i + 1] = x_1[i] + dt * v_1[i];
            v_1[i + 1] = v_1[i] + dt * (-object.w * object.w) * x_1[i];
        }
        for (int i = 0; i < n - 1; i++) {
            object.x[i + 1] = object.x[i] + dt/2 * (object.v[i] + v_1[i + 1]);
            object.v[i + 1] = object.v[i] + dt/2 * (-object.w * object.w)*(object.x[i] + x_1[i + 1]);
        }
    }

    HoynsScheme& operator=(HoynsScheme const &src) = delete;

    HoynsScheme(HoynsScheme const &src) = delete;
};

class EulersMethod{
public:
    EulersMethod() {}

    void count(Object &object){  
        for (int i = 0; i < n - 1; i++){
            object.x[i + 1] = object.x[i] + dt * object.v [i];
            object.v[i + 1] = object.v[i] + dt * (-object.w*object.w) * object.x[i]; 
        }
    }

    EulersMethod& operator=(EulersMethod const &src) = delete;

    EulersMethod(EulersMethod const &src) = delete;

    ~EulersMethod() {}
};

class EulersMethodKahan{
public:
    EulersMethodKahan() {}

    void count(Object &object){
        double c_x = 0;
        double c_v = 0;
        for (int i = 0; i < n - 1; i++){
            object.x[i + 1] = object.x[i] + dt * object.v [i] - c_x;
            c_x = (object.x[i + 1] - object.x[i]) - (dt * object.v [i] - c_x);
            object.v[i + 1] = object.v[i] + dt * (-object.w*object.w) * object.x[i] - c_v;
            c_v = (object.v[i + 1] - object.v[i]) - (dt * (-object.w*object.w) * object.x[i] - c_v);
        }
    }

    EulersMethodKahan& operator=(EulersMethod const &src) = delete;

    EulersMethodKahan(EulersMethodKahan const &src) = delete;

    ~EulersMethodKahan() {}
};

class FileOutput{
public:
    void write(Object &object){
        ofstream out;
        out.open("/home/starman/CLionProjects/EulersMethod/preciseSolution.txt");
        if (out.is_open())
        {
            for (int i = 0; i < n; i++){
                out <<  dt * i  << ' ' << object.x_0*cos(object.w*dt*i) << ' ' <<  object.x_0 * sin(object.w*dt*i) * object.x_0 *sin(object.w*dt*i) / 2 + object.w * object.w *
                object.x_0 * cos(object.w*dt*i) * object.x_0 * cos(object.w*dt*i) / 2 << ' ' << -sin(object.w*dt*i) <<'\n';
            }
        }
        out.close();

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

class Tests{
public:
    void SimpleMethodTest(){
        Object object;
        HoynsSchemeKahan hoynsSchemeKahan;
        hoynsSchemeKahan.count(object);
        object.measureEnergy();
        FileOutput file;
        file.write(object);
    }
    void KahanError() {
        ofstream out;
        out.open("/home/starman/CLionProjects/EulersMethod/TimeError.txt");
        double dt_loc = dt;
        int nn = 500;
        for (int i = 0; i < nn - 2; i++) {
            double id = i;
            dt_loc = dt - dt * id / nn;
            n = floor(duration / dt_loc);
            Object object1, object2, object3;
            HoynsScheme hoynsScheme;
            HoynsSchemeKahan hoynsSchemeKahan;
            hoynsSchemeKahan.count(object1);
            hoynsScheme.count(object2);
            double nd = n;
            for (int i = 0; i < n; i++) {
                object3.x[i] = abs(object1.x[i] - object2.x[i]);
                //object3.v[i] = object1.v[i] - object2.v[i];
                //object3.energy[i] = object1.energy[i] - object2.energy[i];
            }    
            if (out.is_open()) {
                out << dt_loc << ' ' << KahanSum(object3.x, n) / nd << '\n';
            }
        }
        out.close();
    }
        void TimeReverse(){
            Object object, object2;

            HoynsSchemeKahan hoynsSchemeKahan;
            hoynsSchemeKahan.count(object);
            object.measureEnergy();
            dt = -dt;
            object2.x[0] = object.x[n-1];
            object2.v[0] = object.v[n-1];
            hoynsSchemeKahan.count(object2);
            object2.measureEnergy();


            ofstream out;
            out.open("/home/starman/CLionProjects/EulersMEthod/timeReverse.txt");
            if (out.is_open())
            {
                for (int i = 0; i < n; i++){
                    out <<  -dt * (n-i-1)  << ' ' << object2.x[i]<< ' ' << object2.energy[i] << ' ' << object2.v[i] <<'\n';
                }
            }
            out.close();
            dt = abs(dt);
        }
};

int main() {
    Tests test;
    //test.SimpleMethodTest();
    //test.TimeReverse();
    test.KahanError();
    return 0;
}
