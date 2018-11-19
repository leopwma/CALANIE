#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

using namespace std;

/** Prototypes ******************************************************************/

void make_Cijkl(double Cijkl[3][3][3][3], double Cij[6][6]);
void make_Sijkl(double Sijkl[3][3][3][3], double Sij[6][6]);

void read_input(double Cij[6][6], double Sij[6][6], double Pij[3][3]);
/********************************************************************************/

int main(int argc, char **argv){

    double Cij[6][6];
    double Sij[6][6];
    double P[3][3]; //dipole tensor

    read_input(Cij, Sij, P);

    double Cijkl[3][3][3][3];
    double Sijkl[3][3][3][3];
    make_Cijkl(Cijkl, Cij); 
    make_Sijkl(Sijkl, Sij); 

    cout << setprecision(16) ;

    cout << "P11 = " << P[0][0] << " eV \n";
    cout << "P22 = " << P[1][1] << " eV \n";
    cout << "P33 = " << P[2][2] << " eV \n";
    cout << "P12 = " << P[0][1] << " eV \n";
    cout << "P23 = " << P[1][2] << " eV \n";
    cout << "P31 = " << P[2][0] << " eV \n";
    cout << "\n";

    double Omega[3][3]; //relaxation volume tensor
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            Omega[i][j] = 0e0;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    Omega[i][j] += Sijkl[i][j][k][l]*P[k][l];

    cout << "Omega11 = " << Omega[0][0] << " A^3 \n";
    cout << "Omega22 = " << Omega[1][1] << " A^3 \n";
    cout << "Omega33 = " << Omega[2][2] << " A^3 \n";
    cout << "Omega12 = " << Omega[0][1] << " A^3 \n";
    cout << "Omega23 = " << Omega[1][2] << " A^3 \n";
    cout << "Omega31 = " << Omega[2][0] << " A^3 \n";
    cout << "\n";

}            
               
void make_Cijkl(double C[3][3][3][3], double C_Voigt[6][6]){

    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            for (int k = 0; k < 3; ++k){
                for (int l = 0; l < 3; ++l){
                    C[i][j][k][l] = 0e0;
                }
            }
        }
    }
    int p, q, r, s;
    for (int i = 0; i < 6; ++i){
        if (i == 0){
            p = 0;
            q = 0;
        } else if (i == 1){
            p = 1;
            q = 1;
        } else if (i == 2){
            p = 2;
            q = 2;
        } else if (i == 3){
            p = 1;
            q = 2;
        } else if (i == 4){
            p = 2;
            q = 0;
        } else if (i == 5){
            p = 0;
            q = 1;
        }
        for (int j = 0; j < 6; ++j){
            if (j == 0){
                r = 0;
                s = 0;
            } else if (j == 1){
                r = 1;
                s = 1;
            } else if (j == 2){
                r = 2;
                s = 2;
            } else if (j == 3){
                r = 1;
                s = 2;
            } else if (j == 4){
                r = 2;
                s = 0;
            } else if (j == 5){
                r = 0;
                s = 1;
            }
            if (i < 3 && j < 3){
                C[p][q][r][s] = C_Voigt[i][j];
                C[q][p][r][s] = C_Voigt[i][j];
                C[p][q][s][r] = C_Voigt[i][j];
                C[q][p][s][r] = C_Voigt[i][j];
            } else if ((i < 3 && j < 6) || (i < 6 && j < 3)){
                C[p][q][r][s] = C_Voigt[i][j];
                C[q][p][r][s] = C_Voigt[i][j];
                C[p][q][s][r] = C_Voigt[i][j];
                C[q][p][s][r] = C_Voigt[i][j];
            } else {
                C[p][q][r][s] = C_Voigt[i][j];
                C[q][p][r][s] = C_Voigt[i][j];
                C[p][q][s][r] = C_Voigt[i][j];
                C[q][p][s][r] = C_Voigt[i][j];
            }
        }
    }
}

void make_Sijkl(double S[3][3][3][3], double S_Voigt[6][6]){

    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            for (int k = 0; k < 3; ++k){
                for (int l = 0; l < 3; ++l){
                    S[i][j][k][l] = 0e0;
                }
            }
        }
    }
    int p, q, r, s;
    for (int i = 0; i < 6; ++i){
        if (i == 0){
            p = 0;
            q = 0;
        } else if (i == 1){
            p = 1;
            q = 1;
        } else if (i == 2){
            p = 2;
            q = 2;
        } else if (i == 3){
            p = 1;
            q = 2;
        } else if (i == 4){
            p = 2;
            q = 0;
        } else if (i == 5){
            p = 0;
            q = 1;
        }
        for (int j = 0; j < 6; ++j){
            if (j == 0){
                r = 0;
                s = 0;
            } else if (j == 1){
                r = 1;
                s = 1;
            } else if (j == 2){
                r = 2;
                s = 2;
            } else if (j == 3){
                r = 1;
                s = 2;
            } else if (j == 4){
                r = 2;
                s = 0;
            } else if (j == 5){
                r = 0;
                s = 1;
            }
            if (i < 3 && j < 3){
                S[p][q][r][s] = S_Voigt[i][j];
                S[q][p][r][s] = S_Voigt[i][j];
                S[p][q][s][r] = S_Voigt[i][j];
                S[q][p][s][r] = S_Voigt[i][j];
            } else if ((i < 3 && j < 6) || (i < 6 && j < 3)){
                S[p][q][r][s] = S_Voigt[i][j]/2e0;
                S[q][p][r][s] = S_Voigt[i][j]/2e0;
                S[p][q][s][r] = S_Voigt[i][j]/2e0;
                S[q][p][s][r] = S_Voigt[i][j]/2e0;
            } else {
                S[p][q][r][s] = S_Voigt[i][j]/4e0;
                S[q][p][r][s] = S_Voigt[i][j]/4e0;
                S[p][q][s][r] = S_Voigt[i][j]/4e0;
                S[q][p][s][r] = S_Voigt[i][j]/4e0;
            }
        }
    }
}

void read_input(double Cij[6][6], double Sij[6][6], double Pij[3][3])
{
    ifstream in_file_1("input_Pij");
    ifstream in_file_2("input_elastic");


    if (in_file_1) {
        cout << "Start reading input_data ... \n\n";
    } else {
        cout << "ERROR: You need to have an input_Pij file \n\n";
        exit(1);
    }


    while (!in_file_1.eof()){

        char variable[256];
        memset(variable, 0, 256);

        char value[256];
        memset(value, 0, 256);

        in_file_1 >> variable >> value;
        if (strcmp(variable, "P11") == 0) Pij[0][0] = atof(value);
        if (strcmp(variable, "P22") == 0) Pij[1][1] = atof(value);
        if (strcmp(variable, "P33") == 0) Pij[2][2] = atof(value);
        if (strcmp(variable, "P12") == 0) Pij[0][1] = atof(value);
        if (strcmp(variable, "P23") == 0) Pij[1][2] = atof(value);
        if (strcmp(variable, "P31") == 0) Pij[2][0] = atof(value);
    }
    
    Pij[1][0] = Pij[0][1];
    Pij[2][1] = Pij[1][2];
    Pij[0][2] = Pij[2][0];

    in_file_1.close();

    if (in_file_2) {
        cout << "Start reading input_elastic ... \n\n";
    } else {
        cout << "ERROR: You need to have an input_elastic file \n\n";
        exit(1);
    }
   
    string line;  //dummy
    getline(in_file_2,line); //read comment line
    getline(in_file_2,line); //read comment line

    for (int i = 0; i < 6; ++i)
         in_file_2 >> Cij[i][0] >> Cij[i][1] >> Cij[i][2] >> Cij[i][3] >> Cij[i][4] >> Cij[i][5];
    for (int i = 0; i < 6; ++i)
         in_file_2 >> Sij[i][0] >> Sij[i][1] >> Sij[i][2] >> Sij[i][3] >> Sij[i][4] >> Sij[i][5];

    
    double eVperAcubic_to_GPa = 160.21766208;

    for (int i = 0; i < 6; ++i){
        for (int j = 0; j < 6; ++j){
            Cij[i][j] /= eVperAcubic_to_GPa;
            Sij[i][j] *= eVperAcubic_to_GPa;
        }
    } 

    in_file_2.close();

}

