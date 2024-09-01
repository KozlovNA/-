#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <memory>

using namespace std;
#include "json.hpp"
#include "json_fwd.hpp"
using json = nlohmann::json;

// parameters of simulation
double duration;
double dt;
int n;
string OUTPATH;
void init(char* config){
    std::ifstream i(config, std::ifstream::binary);
        json j;
        i >> j;
    duration = j["StartCondition"]["duration"];
    dt = j["StartCondition"]["dt"];
    n = floor(duration/dt);
    OUTPATH = j["OUTPATH"];
}


//-----------MATERIAL DOT OBJECT-----------//
//--contains: start conditions (x_0, v_0)--//
//------------angular frequency (w)--------//
//------------x,v,energy arrays------------//
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
            std::ifstream i("config.json");
            json j;
            i >> j;
            *x = j["StartCondition"]["x_0"];
            *v = j["StartCondition"]["v_0"];
            w = j["StartCondition"]["w"];
        }

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


//---Kahan Summation for array of arguments---//
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


//----------------Kahan Summation for 2 arguments----------------//
//---needs to initiate double c variable for error remembering---//
double KahanSum(double input, double term, double c) {
        double y = term - c;
        double t = input + y;
        c = (t - input) - y;
    return t;
}


//---------------------------------------//
//------NUMERIC CALCULATION METHODS------//
//---------------------------------------//
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

//---------GENERIC--ALGORYTHMS--------//
class PhysPendState{
public:
    array<double, 2> state;
    double g;
    double l;
    PhysPendState(){
        std::ifstream i("config.json", std::ifstream::binary);
        json j;
        i >> j;
        state[0] = j["PhysicalPendulum"]["teta_0"];
        state[1] = j["PhysicalPendulum"]["d(teta_0)/dt"];
        g = j["PhysicalPendulum"]["g"];
        l = j["PhysicalPendulum"]["lambda"];
    }
    PhysPendState(double x, double v){
        PhysPendState st;
        g = st.g;
        l = st.l;
        state[0] = x;
        state[1] = v;
    }
    PhysPendState(PhysPendState const &single_state){
        PhysPendState st;
        g = st.g;
        l = st.l;
        state[0] = single_state.state[0];
        state[1] = single_state.state[1];
    }

    array<double, 2> f(PhysPendState &state){
        array<double, 2> result;
        result[0] = state.state[1];
        result[1] = -g/l*sin(state.state[0]);
        return result;
    }

    array<double, 2> sum(array<double, 2> first, array<double, 2> second){
        array<double, 2> result;
        result[0] = first[0]+second[0];
        result[1] = first[1]+second[1];
        return result;
    }

    array<double, 2> prod(array<double, 2> first, double second){
        array<double, 2> result;
        result[0] = first[0]*second;
        result[1] = first[1]*second;
        return result;
    }
};


template<typename State>
class GenericEuler{
public:
    State step(State &state){
        State new_state;
        new_state.state = state.sum(state.state, state.prod(state.f(state), dt));
        return new_state;
    } 
};


template<typename State, typename Method>
void solver(){
    Method method;
    vector<State> states;
    State initstate;
    states.push_back(initstate);
    for (int i = 0; i < n; i++){
        State new_state(method.step(states.back())); 
        states.push_back(new_state);
    }
    ofstream out;
        out.open(OUTPATH);
        if (out.is_open())
        {
            for (int i = 0; i < n; i++){
                out <<  dt * i  << ' ' << states[i].state[0]<< ' ' << states[i].state[1] <<'\n';
            }
        }
        out.close();
}

//--------------------------------------//
class RungeKutta{
public:
    RungeKutta() {}
    ~RungeKutta(){}

    double *f(Object &object, int i){
        static double res[2] = {0, 0};
        res[0] = object.v[i];
        res[1] = -object.x[i]*object.w*object.w;
        return res;
    }

    double *k_1(Object &object, int i){
        static double res[2] = {0, 0};
        res[1] = f(object, i)[1]*dt;
        res[0] = f(object, i)[0]*dt;
        return res;
    }
    double *k_2(Object &object, int i){
        static double res[2] = {0, 0};
        Object object2;
        object2.x[i] = object.x[i] + 1.0/2.0 * k_1(object, i)[0];
        object2.v[i] = object.v[i] + 1.0/2.0 * k_1(object, i)[1];
        res[0] = f(object2, i)[0]*dt;
        res[1] = f(object2, i)[1]*dt;
        return res;
    }
    double *k_3(Object &object, int i){
        static double res[2] = {0, 0};
        Object object2;
        object2.x[i] = object.x[i] + 1.0/2.0 * k_2(object, i)[0];
        object2.v[i] = object.v[i] + 1.0/2.0 * k_2(object, i)[1];
        res[0] = f(object2, i)[0]*dt;
        res[1] = f(object2, i)[1]*dt;
        return res;
    }
    double *k_4(Object &object, int i){
        static double res[2] = {0, 0};
        Object object2;
        object2.x[i] = object.x[i] + k_3(object, i)[0];
        object2.v[i] = object.v[i] + k_3(object, i)[1];
        res[0] = f(object2, i)[0]*dt;
        res[1] = f(object2, i)[1]*dt;
        return res;
    }

    void count(Object &object){  
        for (int i = 0; i < n; i++){
            object.x[i+1] = object.x[i] + 1.0/6.0 * (k_1(object, i)[0] + 2*k_2(object, i)[0] + 2*k_3(object, i)[0] + k_4(object, i)[0]); 
            object.v[i+1] = object.v[i] + 1.0/6.0 * (k_1(object, i)[1] + 2*k_2(object, i)[1] + 2*k_3(object, i)[1] + k_4(object, i)[1]);
        }
    }

    RungeKutta& operator=(RungeKutta const &src) = delete;

    RungeKutta(RungeKutta const &src) = delete;
};


//---------------------------------------//
//-------SCRIPTS FOR DATA OUTPUT---------//
//---------------------------------------//
class FileOutput{
public:
    void write(Object &object){
        ofstream out;
        out.open("/home/starman/CLionProjects/RangiCut/preciseSolution.txt");
        if (out.is_open())
        {
            for (int i = 0; i < n; i++){
                out <<  dt * i  << ' ' << object.x_0*cos(object.w*dt*i) << ' ' <<  object.x_0 * sin(object.w*dt*i) * object.x_0 *sin(object.w*dt*i) / 2 + object.w * object.w *
                object.x_0 * cos(object.w*dt*i) * object.x_0 * cos(object.w*dt*i) / 2 << ' ' << -sin(object.w*dt*i) <<'\n';
            }
        }
        out.close();

        out.open("/home/starman/CLionProjects/RangiCut/RungeKuttData.txt");
        if (out.is_open())
        {
            for (int i = 0; i < n; i++){
                    out <<  dt * i  << ' ' << object.x[i]<< ' ' << object.energy[i] << ' ' << object.v[i] <<'\n';
            }
        }
        out.close();
    }
};

//-------VARIOUS BASIC TESTS AND EXPERIMANTS---------//
//--Simple Method Test------------calculate trajectory, using object's start conditions
//--KahanError--------------------estimates errror that occures depending on machine epsilon
//--Time Reverse------------------calculates trajectory in strait direction and then backwards
//--Time Reverse Error------------estimates error depending on scale of dt 
class Tests{
public:


    void SimpleMethodTest(){
        Object object;
        RungeKutta method;
        method.count(object);
        object.measureEnergy();
        FileOutput file;
        file.write(object);
    }


    void KahanError() {
        ofstream out;
        out.open("/home/starman/CLionProjects/RangiCut/TimeError.txt");
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
            Tests test;
            test.SimpleMethodTest();
            Object object, object2;

            HoynsSchemeKahan method;
            method.count(object);
            object.measureEnergy();
            dt = -dt;
            object2.x[0] = object.x[n-1];
            object2.v[0] = object.v[n-1];
            method.count(object2);
            object2.measureEnergy();


            ofstream out;
            out.open("/home/starman/CLionProjects/RangiCut/timeReverse.txt");
            if (out.is_open())
            {
                for (int i = 0; i < n; i++){
                    out <<  -dt * (n-i-1)  << ' ' << object2.x[i]<< ' ' << object2.energy[i] << ' ' << object2.v[i] <<'\n';
                }
            }
            out.close();
            dt = abs(dt);
        }

        
        void TimeReverseError(){
        ofstream out;
        out.open("/home/starman/CLionProjects/RangiCut/TimeReverseError.txt");
        double dt_loc = dt;
        int nn = 500;
        for (int i = 0; i < nn; i++) {
            double id = i;
            dt_loc = dt - dt * id / nn;
            n = floor(duration / dt_loc);
            Object object1, object2, object3;
            EulersMethod method;
            method.count(object1);
            object1.measureEnergy();
            dt_loc = -dt_loc;
            object2.x[0] = object1.x[n];
            object2.v[0] = object1.v[n];
            method.count(object2);
            object2.measureEnergy();
            dt_loc = -dt_loc;
            for (int j = 0; j < n; j++) {
                object3.x[j] = abs(object1.x[j] - object2.x[n-j]);
                //object3.v[i] = object1.v[i] - object2.v[i];
                object3.energy[j] = abs(object1.energy[j] - object2.energy[n-j]);
            } 
            double nd = n;   
            if (out.is_open()) {
                out << dt_loc << ' ' << KahanSum(object3.x, n) / nd <<  ' ' << KahanSum(object3.energy, n) / nd << '\n';
            }
        }
        out.close();
        n = floor(duration/dt);
        }
};

int main(int argc, char* argv[]) {
    init("config.json");
    //Tests test;
    //test.SimpleMethodTest();
    //test.TimeReverse();
    //test.KahanError();
    //test.TimeReverseError();
    solver<PhysPendState, GenericEuler<PhysPendState>>();
    return 0;
}