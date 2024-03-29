/*The program below consider the problems of particle fling over the semi-infinite potential of rectangular barriers
 * located at a distance from each other much more than their characteristic size (a) at different initial values of
 * energy (E). Reflection consider only in cases of returning particle after the first barrier.
 */

#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <string>
#include <omp.h>

const double eV = 1.602176634e-19; //International System of Units.
const double m = 1.6726219e-27; //Mass of particle.
const double h = 1.054571817e-34; //Planck's constant.
const double U = 19*eV; //Default height of the barrier.
const double eps = 0.1; //Accuracy.

//Function counts transmission coefficient (D).
double Passing(double k1, double k2) {
    return 4 * k1 * k2 / std::pow(k1 + k2, 2);
}

//Function computes the total reflectance of particles, which returns after first barrier.
bool Reflection(double E, double& Refl, int l, int sign) {
    double D, R = 1.0;
    while (R > eps && (E > U || l > 0) && !std::isnan(E)) {
        double k1 = std::sqrt(2.0 * m * E) / h;
        double k2 = std::sqrt(2.0 * m * (E - U)) / h;
        D = Passing(k1, k2);
        E *= D;
        R = 1 - D;
        l += sign;
        if (Reflection(E, R, l, -sign) && !std::isnan(R))
            Refl += R;
    }
    return (l == 0);
}

/*Function fills vector of pairs (E, R) in order in parallel.
 * May be it would be better to fill the std::vector in no particular order
 * then use std::sort() with comparator, but it works - stay out.*/
std::vector<std::pair <double, double>> DataSetCreation(std::vector<double> E) {
    std::vector<std::pair<double, double>> EnRef;
    #pragma omp parallel
    {
        std::vector<std::pair<double, double>> EnRef_private;
        #pragma omp for nowait schedule(static)
        for (unsigned i = 0; i < E.size(); i++) {

            //Initials.
            double R = 0;
            int l = 0;
            int direction = 1;

            double Ebuf = E[i] * eV;
            Reflection(Ebuf, R, l, direction); //tmp -- just a variable for calling a function;
            EnRef_private.emplace_back(std::move(std::make_pair(E[i], R)));
        }
        #pragma omp for schedule(static) ordered
        for (unsigned i = 0; i < omp_get_num_threads(); i++) {
            #pragma omp ordered
            EnRef.insert(EnRef.end(), EnRef_private.begin(), EnRef_private.end());
        }
    }
    return EnRef;
}

void DataFileCreation (std::string& DataType, std::vector<std::pair<double, double>>& xx){
    //For reading created files via Matlab use command: M = dlmread('/PATH/file'); xi = M(:,i);
    std::ofstream fout;
    fout.open(DataType);
    for(unsigned i = 0; i < xx.size(); i++)
        fout << xx[i].first << '\t' << xx[i].second << std::endl;
    fout.close();
}

//Function creates plots via GNUPlot.
void plot(std::string& name, std::string& data, std::string xlabel, std::string ylabel, std::string title) {
    //if you have problems with ".svg", you have to change ".svg" to ".pdf" in strings bellow.
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp)
        throw std::runtime_error("Error opening pipe to GNU plot.");
    else {
        std::vector<std::string> stuff = {"set term svg",
                                          "set out \'" + name + ".svg\'",
                                          "set xlabel \'" + xlabel + "\'",
                                          "set ylabel \'" + ylabel + "\'",
                                          "set grid xtics ytics",
                                          "set title \'" + title + "\'",
                                          "plot \'" + data + "\'using 1:2 with lines lw 2 lt rgb 'blue',\
                                          \'" + data + "\' using 1:2 lw 1 lt rgb 'orange' ti \'Nodes\'",
                                          "set key box top right",
                                          "set terminal wxt",
                                          "set output",
                                          "replot"};
        for (const auto &it : stuff)
            fprintf(gp, "%s\n", it.c_str());
        pclose(gp);
    }
}

int main() {
    double FirstTask = 20;
    double LastTask = 200;
    std::vector<double> E(LastTask - (FirstTask - 1)); //Vector of problems.
    std::generate(E.begin(), E.end(), [&] {return FirstTask++;}); //Range of potential barrier height.
    std::vector<std::pair<double, double>> data = std::move(DataSetCreation(E));
    std::string DataType = "Reflection";
    DataFileCreation(DataType, data);
    plot(DataType, DataType, "E, eV", "R = R(E)", DataType);
    return 0;
}
