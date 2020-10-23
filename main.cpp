/*The program below consider the problems of particle fling over the semi-infinite potential of rectangular barriers
 * located at a distance from each other much more than their characteristic size (a) at different initial values of
 * energy (E).
 */

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <string>
#include <tuple>


const double Ev = 1.602176634e-19;
const double m = 1.6726219e-27; //Mass of particle.
const double h = 1.054571817e-34; //Planck's constant.
const double a = 1e-10; //characteristic size.
const double U = 19*Ev; //Default height of the barrier.


void pass(std::complex<double>& jPssd, double& D, double k1, double k2) {
    //Function counts the passed flow (jPssd) and transmission coefficient (D) for every task.
    const auto i = std::complex<double>(0, 1);
    std::complex<double> d3pow2 = -(4.0*std::pow(k2, 2) * std::exp(-a * k1 * (2.0*i)) * std::exp(a*k2*(2.0*i)) *
                                    std::pow((k1 * (2.0*i) - k2 * i), 2)) / std::pow((k1*k2 - std::pow(k1, 2) *
                                    std::exp(a * k2 * (2.0*i)) - std::pow(k2, 2) * std::exp(a * k2 * (2.0*i)) +
                                    std::pow(k1, 2) + 2.0 * k1 * k2 * std::exp(a * k2 * (2.0*i))), 2);
    std::complex<double> jInc = h * k1 / m; //Incoming flow for current barrier.
    jPssd = h * k1 / m / d3pow2;
    D = std::abs(jPssd / jInc);
}

std::tuple<double, double, unsigned long> Reflection(double E) {
    double jInc, D;
    std::complex<double> jPssd;
    bool flag = true;
    unsigned long l = 0; //Number of current barrier;
    double Ebuf = E;
    E *= Ev;
    do {
        double k1 = std::sqrt(2.0 * m * E) / h;
        double k2 = std::sqrt(2.0 * m * (E - U)) / h;
        pass(jPssd, D, k1, k2);
        //buf = E;
        E *= D;
        if (flag == 1) {
            jInc = h * k1 / m; //Incoming flow.
            flag = false;
        }
        l++;
    } while (E > U); //The case of tunneling is not considered.
    double R = 1 - std::abs(jPssd / jInc);
    return std::make_tuple(Ebuf, R, l);
}

void dataset (std::string& DataType, std::vector<std::tuple<double, double, unsigned long>>& xx){
    //Function creates data-files for plots with key.
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
    for (const auto& it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}

int main() {
    double e0 = 20; //Energy of the particle in the first case.
    double eF = 200; //Energy (Ev) of the particle in the last case.
    std::vector<double> E(eF - (e0 - 1)); //Number of problems.
    std::generate(E.begin(), E.end(), [&] {return e0++;}); //Scatter of potential barrier height.
    std::vector<std::tuple<double, double, unsigned long>> EnRef;
    for (unsigned i = 0; i < E.size(); i++)
        EnRef.push_back(Reflection(E[i]));
    std::string DataType = "Reflection";
    dataset(DataType, EnRef);
    plot("Reflection", "Reflection", "Energy, Ev", "Reflection", 1, "R = R(E)");
    plot("Barriers", "Reflection", "Energy, Ev", "Number of potential barriers", 2, "Barriers passed");
    return 0;
}
