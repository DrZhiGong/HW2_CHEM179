//Dear Xiao,This assignment was a lot of work in programming and using armadillo, I used Visual Studio as the compiler and followed the settings in the reference link attached at the end, but there were still some problems, I had some path conflicts on my own laptop, so I made sure that the code worked as best as I could, and luckily the output was the same as yours, but with a different decimal place. is the same as yours but with different retained decimal places, more importantly I refer to ChatGPT's output format for this codeframe and refer to some articles and tweets from Extended Huckel Theory, the link will be attached as well, in addition you can let me know if you run into any problems when testing the code, the score does matter to me, especially since I'm This score is really important to me, especially when I apply for PhD position at the end of the year. I hope all is well in your life on the research path! I have added you to Github. Reference website: 1>https://www.researchgate.net/publication/343193259_AN_IMPROVED_FAMILY_OF_BLOCK_METHODS_BASE_ON_THE_EXTENDED_TRAPEZOIDAL_RULE_OF_SECOND_KIND_AND_THEIR_APPLICATIONS_AN_IMPROVED_FAMILY_OF_BLOCK_METHODS_BASE_ON_THE_EXTENDED_TRAPEZOIDAL_RULE_OF_SECOND_KIND 2>https://math.stackexchange.com/questions/2891298/derivation-of-2d-trapezoid-rule 3>https://blog.csdn.net/qq_27278957/article/details/78075632 4>https://blog.csdn.net/dkink/article/details/4787686 5>https://blog.csdn.net/zhebushibiaoshifu/article/details/127123511?ops_request_misc=&request_id=&biz_id=102&utm_term=%E5%A6%82%E4%BD%95%E4%BD%BF%E7%94%A8armadillo&utm_medium=distribute.pc_search_result.none-task-blog-2~all~sobaiduweb~default-2-127123511.142^v99^pc_search_result_base3&spm=1018.2226.3001.4187 6>https://youtu.be/YsRe2PMCbIw?si=7Jp_kg5boB-0qK5l 7>ChatGPT
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

double M_PI = 3.1415926;// ATTENTION: The number of Pi is a approximation, may be not fine enough.

// For a simple clarification about the single atom Gaussian Function.
double gaussiansss(double x, double Xcenter, double alpha, double l) {
    return pow(x - Xcenter, l) * exp(-alpha * pow(x - Xcenter, 2));
}
// n!
int n1(int n) {
    int res = 1, i;
    for (i = 1; i <= n; i++) {
        res *= i;
    }
    return res;
}

// n!!
int n11(int n) {
    int i; double res = 1.0;
    for (i = n; i >= 1; i -= 2) {
        res *= i;
    }
    return res;
}

// BioCoefficient Finding.
int BioCoefficient(int m, int n) {
    return n1(m) / (n1(n) * n1(m - n));
}
// Single dimension integration as CHEM179 LAB2 shows.
double D1interg(double x, double Xcenter, double Ycenter, double alpha, double beta, double X1, double Y2) {
    return gaussiansss(x, Xcenter, alpha, X1) * gaussiansss(x, Ycenter, beta, Y2);
}

// For xyz system.
double XYZcombin(double Xcenter, double Ycenter, double alpha, double beta, double x1, double y1) {

    double p = exp(-alpha * beta * pow(Xcenter - Ycenter, 2) / (alpha + beta));
    double Xp = (alpha * Xcenter + beta * Ycenter) / (alpha + beta);
    double disum = 0;

    for (int i = 0; i < x1 + 1; i++) {
        double insideones = 0;
        for (int j = 0; j < y1 + 1; j++) {
            if ((i + j) % 2 == 0) {
                double summand = BioCoefficient(x1, i) * BioCoefficient(y1, j) * n11(i + j - 1) * pow(Xp - Xcenter, x1 - i) * pow(Xp - Ycenter, y1 - j) / pow(2 * (alpha + beta), (i + j) / 2);
                insideones += summand;
            }
        }
        disum += insideones;
    }
    return p * sqrt(M_PI / (alpha + beta)) * disum;
}

// Extended Trapezoidal Rule
double ETRule(double a, double b, int n, double Xcenter, double Ycenter, double alpha, double beta, double x1, double y1)
{
    double descent = (b - a) / (n - 1);
    double value = 0.5 * (D1interg(a, Xcenter, Ycenter, alpha, beta, x1, y1) + D1interg(b, Xcenter, Ycenter, alpha, beta, x1, y1));

    for (int k = 2; k < n; k++) {
        value += D1interg(a + descent * (k - 1), Xcenter, Ycenter, alpha, beta, x1, y1);
    }
    value *= descent;
    return value;
}

// Extended Trapezoidal Rule 2 Ref. website:https://www.researchgate.net/publication/343193259_AN_IMPROVED_FAMILY_OF_BLOCK_METHODS_BASE_ON_THE_EXTENDED_TRAPEZOIDAL_RULE_OF_SECOND_KIND_AND_THEIR_APPLICATIONS_AN_IMPROVED_FAMILY_OF_BLOCK_METHODS_BASE_ON_THE_EXTENDED_TRAPEZOIDAL_RULE_OF_SECOND_KIND
double ETRule2(double a, double b, double n, double Xcenter, double Ycenter, double alpha, double beta, double x1, double y1) {
    double x, m=0, summary, delta;
    static float s;
    int it, j;
    if (n == 1) {
        return (s = 0.5 * (b - a) * (D1interg(a, Xcenter, Ycenter, alpha, beta, x1, y1) + D1interg(b, Xcenter, Ycenter, alpha, beta, x1, y1)));
    }
    else {
        for (it = 1, j = 1; j < n - 1; j++) {
            it <<= 1;
            delta = (b - a) / m;
            x = a + 0.5 * delta;
            m = it;
            for (summary = 0.0, j = 1; j <= it; j++, x += delta) {
                summary += D1interg(x, Xcenter, Ycenter, alpha, beta, x1, y1);
            }
            s = 0.5 * (s + (b - a) * summary / m);
        }

        return s;
    }
}

// Two Atom calculations.
double D1intergoverlap(double Xcenter, double Ycenter, double alpha, double beta, double x1, double y1) {
    double a, b, stdA, stdB;
    stdA = sqrt(1 / (2 * alpha));
    stdB = sqrt(1 / (2 * beta));
    if ((Xcenter - 4 * stdA) <= (Ycenter - 4 * stdB)) {
        a = Xcenter - 4 * stdA;
    }
    else {
        a = Ycenter - 4 * stdB;
    }
    if ((Xcenter + 4 * stdA) >= (Ycenter + 4 * stdB)) {
        b = Xcenter + 4 * stdA;
    }
    else {
        b = Ycenter + 4 * stdB;
    }
    double n = 10000; 
    return ETRule(a, b, n, Xcenter, Ycenter, alpha, beta, x1, y1);
}

// Just like HW1`s file reading part.
int Counter(string file_name) {
    ifstream inputfile;
    inputfile.open(file_name);

    string line;
    int atomnum = 0;

    while (getline(inputfile, line)) {
        stringstream lineStream(line);
        while (getline(lineStream, line, ' ')) {
            atomnum++;
        }
    }

    inputfile.close();
    return atomnum;
}

//  Identify all possible combinations of angular momentum components for a given atomic shell.
vector<vector<double>> findTriplets(int L) {
    vec bank = linspace<vec>(L, 0, L + 1);
    vector<vector<double>> triplets;
    int n = bank.n_elem;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                if (bank[i] + bank[j] + bank[k] == L) {
                    vector<double> triplet = { bank(i),bank(j),bank(k) };
                    triplets.push_back(triplet);
                }
            }
        }
    }
    return triplets;
}



int main() {
    string file_name = "\Myarmadillo\testhw2\sample_input\analytical\3.txt";

    // Counter.
    if (Counter(file_name) < 7) {

        ifstream inputFile(file_name);
        if (!inputFile) {
            cerr << "Error opening file." << endl;
        }

        // Get variables from imput file. Ref.Lab2
        double Xcenter, Ycenter, alpha, beta, x1, y1;
        inputFile >> Xcenter >> alpha >> x1;
        inputFile >> Ycenter >> beta >> y1;


        inputFile.close();  

        // Calculate overlap integral
        double question1res = D1intergoverlap(Xcenter, Ycenter, alpha, beta, x1, y1);
     
        cout << "1d numerical overlap integral between Gaussian functions is " << question1res << endl;
    }//End the Question 1.
    else {
   
        ifstream inputFile(file_name);

        if (!inputFile) {
            cerr << "Error opening file." << endl;
        }

        double Xcenter, Ycenter, Zcenter, Xcenter2, Ycenter2, Zcenter2, alpha, beta, Lbond1, Lbond2;
        inputFile >> Xcenter >> Ycenter >> Zcenter >> alpha >> Lbond1;
        inputFile >> Xcenter2 >> Ycenter2 >> Zcenter2 >> beta >> Lbond2;
        
        inputFile.close(); 


        // Calculate number of subshells/functions for the shell
        double Multiorbital1 = (Lbond1 + 1) * (Lbond1 + 2) / 2;
        double Multiorbital2 = (Lbond2 + 1) * (Lbond2 + 2) / 2;

        cout << "Shell 1 has " << Multiorbital1 << " functions." << endl;
        cout << "This shell info: R(" << Xcenter << ", " << Ycenter << ", " << Zcenter << "), with angular momentum: " << Lbond1 << ", coefficient: " << alpha << endl;
        cout << "Shell 2 has " << Multiorbital2 << " functions." << endl;
        cout << "This shell info: R(" << Xcenter2 << ", " << Ycenter2 << ", " << Zcenter2 << "), with angular momentum: " << Lbond2 << ", coefficient: " << beta << endl;

        vector<vector<double>> Tri1 = findTriplets(Lbond1);
        vector<vector<double>> Tri2 = findTriplets(Lbond2);

        // Overlap Matrix.
        mat OverlapIntegrals(Multiorbital1, Multiorbital2, fill::zeros); 
        for (int i = 0; i < Multiorbital1; i++) { 
            for (int j = 0; j < Multiorbital2; j++) { 

                double lA = Tri1[i][0], mA = Tri1[i][1], nA = Tri1[i][2];
                double lB = Tri2[j][0], mB = Tri2[j][1], nB = Tri2[j][2];

                double atoms3 = 1; 
                atoms3 *= XYZcombin(Xcenter, Xcenter2, alpha, beta, lA, lB); 
                atoms3 *= XYZcombin(Ycenter, Ycenter2, alpha, beta, mA, mB); 
                atoms3 *= XYZcombin(Zcenter, Zcenter2, alpha, beta, nA, nB); 
                OverlapIntegrals(i, j) = atoms3; 
            }
        }
        // Output XYZ 3D matrix.
        OverlapIntegrals.print("Overlap integral between Shell 1 and Shell 2");
   

        cout << "The components of angular momentum (l, m, n) for the matrix column, from top to bottom, are listed sequentially as: (";
        for (int i = 0; i < Tri1.size(); i++) {
            cout << Tri1[i][0] << ", " << Tri1[i][1] << ", " << Tri1[i][2];
            if (i != Tri1.size() - 1) {
                cout << "), (";
            }
            else {
                cout << ")." << endl;
            }
        }
        cout << "The components of angular momentum (l, m, n) for the matrix row, from left to right, are listed sequentially as: (";
        for (int i = 0; i < Tri2.size(); i++) {
            cout << Tri2[i][0] << ", " << Tri2[i][1] << ", " << Tri2[i][2];
            if (i != Tri2.size() - 1) {
                cout << "), (";
            }
            else {
                cout << ")." << endl;
            }
        }//End the Question 2.
    }
    return 0;
}


