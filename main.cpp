#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <string>
#include <tuple>
#include "omp.h"

//Basic constants
const double h = 1.054571817e-34; //Planck's constant divided by 2*pi.
const double e = -1.6021766208e-19; //Electron charge. //?!
const double N_A =  6.02214076e23; //Avogadro constant.
const double a_TF = 0.47e-8; //Bohr radius.

//Incoming atom
const double Z1 = 1; //Number in the periodic table.
const double m_1 = 1.6726219e-24; //Atomic mass.

//Material properties
const double Z2 = 13; //Number in the periodic table.
const double m_a = 26.9815386; //Atomic mass.
const double rho = 2.6989; //Density.
const double N = 1.0 / m_a * N_A * rho; //The number of atoms in one cubic centimeter.
const double m_2 = 1.660539066601e-24 * 26.981539;

//Problem constants
const double m = (m_1*m_2) / (m_1 + m_2);
const double factor = (0.831 * Z1 * Z2 * std::pow(e, 2)) /
                      (2 * std::sqrt(std::pow(Z1, (2.0/3.0)) + std::pow(Z2, (2.0/3.0)))) * a_TF;

//Program constants
const int sample = 1e4;
const double eps = 1.0/sample;

double ScreenedCoulombPotential (double r) {
    return factor / std::pow(r, 2);
}

double EntryCoordinate(double E) {
    return std::sqrt(factor / E);
}

void pass(double& D, double& R, double k1, double k2) {
    D = k1*std::pow((k1 + k2), 2) / (4*std::pow(k1, 2)*k2);
    R += k1*std::pow((k1 - k2), 2) / (4*std::pow(k1, 2)*k2);
}

std::tuple<double, double, unsigned long> Reflection(double E) {
    int l = 0;
    double D, R = 0;
    do {
        double k1 = std::sqrt(2.0 * m * E) / h;
        double a = EntryCoordinate(E);
        double k2 = std::sqrt(2.0 * m * (ScreenedCoulombPotential(a) - E)) / h;
        pass(R, D, k1, k2);
        E *= D;
        }
        l++;
    } while (E > 0); // && std::abs(E - buf) >= eps * buf); //The case of tunneling is not considered.
    double x0 = 0; //Just a start point for grid.
    std::vector<double> grid (sample + 1);
    std::generate(grid.begin(), grid.end(), [&] {return x0 += eps;});
    std::vector<double> spline = CubicSpline(l, Ex, grid);
    return std::make_tuple(Ebuf, R, l);
}

std::vector<double> SimpleSplineMoments(std::vector<double>& f, double step) {
    std::vector<double> moments;
    unsigned n = f.size()-1;
    moments.push_back((3*f[0] - 4*f[1] + f[2])/2.0/step);
    for(unsigned i = 1; i < n; i++)
        moments.push_back((f[i+1] - f[i-1])/2.0/step);
    moments.push_back((-3*f[n] + 4*f[n-1] - f[n-2])/2.0/step);
    return moments;
}

std::vector<double> CubicSpline(std::vector<double>& x, std::vector<double>& f, std::vector<double>& xx){
    //Function returns vector of spline points.
    std::vector<double> yy;
    double step = x[1] - x[0];
    unsigned i = 0;
    std::vector<double> moments = SimpleSplineMoments(f, step);
    for(unsigned j = 0; j < x.size()-1; j++)
        while(xx[i] >= x[j] && xx[i] <= x[j+1]){
            double buf1 = x[j+1] - xx[i];
            double buf2 = xx[i] - x[j];
            yy.push_back(pow(buf1, 2)*(2*buf2+h)*f[j]/pow(h, 3)+
                         pow(buf2, 2)*(2*buf1+h)*f[j+1]/pow(h, 3)+
                         pow(buf1, 2)*buf2*moments[j]/pow(h, 2)-
                         pow(buf2, 2)*buf1*moments[j+1]/pow(h, 2));
            i++;
        }
    return yy;
}

double ConstStepTrapIntegral(std::vector<double>& x, std::vector<double>& f) {
    double I = 0;
    for(unsigned i = 1; i < x.size(); i++)
        I += (f[i-1] + f[i]) / 2.0 * (x[i] - x[i-1]);
    return I;
}



void dataset (std::string& DataType, std::vector<std::tuple<double, double, unsigned long>>& xx){
    //Function creates data-files for plots.
    //For reading created files via Matlab use command: M = dlmread('/PATH/file'); xi = M(:,i);
    std::ofstream fout;
    fout.open(DataType);
    for(unsigned i = 0; i < xx.size(); i++)
        fout << std::get<0>(xx[i]) << '\t' << std::get<1>(xx[i]) << '\t' << std::get<2>(xx[i]) << std::endl;
    fout.close();
}

template <typename T> //The template bellow using for transformation num to str.
std::string toString(T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}

void plot(std::string name, std::string data, std::string xlabel, std::string ylabel, unsigned ycol, std::string title) {
    //Creates plots with data-files and key.
    //Change ".svg" to ".pdf", if you have problems with ".svg".
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp)
        throw std::runtime_error("Error opening pipe to GNU plot.");
    std::vector<std::string> stuff = {"set term svg",
                                      "set out \'" + name + ".svg\'",
                                      "set xlabel \'" + xlabel + "\'",
                                      "set ylabel \'" + ylabel + "\'",
                                      "set grid xtics ytics",
                                      "set title \'" + title + "\'",
                                      "plot \'" + data + "\'using 1:" + toString(ycol+1) + " with lines lw 2 lt rgb 'blue',\
                                      \'" + data + "\' using 1:" + toString(ycol+1) + " lw 1 lt rgb 'orange' ti \'Nodes\'",
                                      "set key box top left",
                                      "set terminal wxt",
                                      "set output",
                                      "replot"};
    for (const auto &it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}

bool comparator(const std::tuple<double, double, unsigned long>& a, const std::tuple<double, double, unsigned long>&  b) {
    //Comparator makes std::sort sort the std::vector of std::pair through the first values in std::pairs.
    return std::get<0>(a) < std::get<0>(b);
}

int main() {
    double e0 = 20; //Energy of the particle in the first case.
    double eF = 200; //Energy (Ev) of the particle in the last case.
    std::vector<double> E(eF - (e0 - 1)); //Number of problems.
    std::generate(E.begin(), E.end(), [&] {return e0++;}); //Scatter of potential barrier height.
    std::vector<std::tuple<double, double, unsigned long>> EnRef;
//#pragma omp parallel for default(none) shared(E, EnRef)
    for (unsigned i = 0; i < E.size(); i++)
        EnRef.push_back(Reflection(E[i]));
    //std::sort(EnRef.begin(), EnRef.end(), comparator); //Sorts vector of energies and reflection coefficients through energies.
    std::string DataType = "Reflection";
    dataset(DataType, EnRef);
    plot("Reflection", "Reflection", "Energy, Ev", "Reflection", 1, "R = R(E)");
    plot("Barriers", "Reflection", "Energy, Ev", "Number of potential barriers", 2, "Barriers passed");
    return 0;
}