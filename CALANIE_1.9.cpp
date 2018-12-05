/********************************************************************************
*
*   Copyright (C) 2018 Culham Centre for Fusion Energy,
*   United Kingdom Atomic Energy Authority, Oxfordshire OX14 3DB, UK
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*
********************************************************************************
*
*   Program: CALANIE
*            CALculation of ANIsotropic Elastic interaction energy of
*            a defect in periodic boundary conditions
*   Version: 1.8
*   Date:    18 Oct 2018
*   Author:  Pui-Wai (Leo) MA
*   Contact: leo.ma@ukaea.uk
*   Address: Culham Centre for Fusion Energy, OX14 3DB, United Kingdom
*
********************************************************************************
* 
*  This program calculate the elastic interaction energy of a defect with its
*  images in periodc boundary conditions. It also output the elastic dipole 
*  tensor and relaxation volume tensor.
*  
*  It requires two input files: input_data & input_elastic.
*
*  "input_data" will be discussed here.
*
*  "input_elastic" contains the information of elastic constant tensor.
*  Format of it is mentioned in "make_input_elastic.py".
*
*  Compilation:
*
*  One can compile the program in 2 ways.
*
*  1) g++ -DORIENTATION CALANIE_version.cpp
*
*     This requires input of 
*     (1) Reference simulation box vectors
*     (2) Lattice unit for reference box vectors
*     (3) Relaxation volumes Omega1 and Omega2
*     (4) Orientation of the defect
*
*     Put (1) to (3) in a file with name "input_data", and also the 
*     "input_elastic". 
*
*     Run the program by
*     ./a.out theta phi
*     where theta and phi are (4) Orientation of the defect
*
*     Detail is in S. L. Dudarev and Pui-Wai Ma, 
*     Physical Review Materials 2 (3), 033602 (2018)
*
*  Or
*
*  2) g++ -DABINITIO -DSTRESSeV CALANIE_version.cpp
*     or
*     g++ -DABINITIO -DSTRESSGPa CALANIE_version.cpp
*
*     This require input of
*     (1) Reference simulation box vectors
*     (2) Lattice units for reference simulation box vectors
*     (3) Total marcoscopic stress of the reference simulation box (in eV or GPa) 
*     (4) Defect simulation box vectors
*     (5) Lattice units for defect simulation box vectors
*     (6) Total marcoscopic stress of the defect simulation box (in eV or GPa)
*
*     Put (1) to (6) in a file with name "input_data"
*
*     Run the program by
*     ./a.out
*
*     Note:
*     
*     -DSTRESSeV means total stresses multiply volume. Unit is in eV. It is 
*     exactly the line of "Total" in vasp "OUTCAR". 
*
*     -DSTRESSGPa means total stresses. Unit is in GPa. It is corresponding to 
*     the line of "kB" in vasp "OUTCAR", where "kB" means kbar = 0.1GPa.
*
********************************************************************************
* Optional:
*    (1) In the correct term, we are doing cubic sum. We can change it into
*        spherical sum by adding -DSPHERICALSUM . 
*    (2) We observed that if we use cubic sum, the correction term tends to
*        zero, if we consider more neighbour cells. One can ignore the 
*        calculation of the correction term by adding -DNOCORRECT. The default 
*        value of Range_Neigh changes from 10 to 30. However, it is much faster.  
*    (3) (Very important) It is strongly recommended to add either 
*        -DSPHERICALSUM or -DNOCORRECT to the compilation, but not both.
********************************************************************************
*
*  A sample "input_data" is given for both cases. "input_data_1" is for 1st 
*  case. "input_data_2" is for 2nd case using -DSTRESSeV. 
*
*  Sign convention of total marcoscopic stress follow VASP. 
*  i.e positive means the simulation cell intends to expand.
*
********************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

using namespace std;

/** Prototypes ******************************************************************/

double delta(int i, int j);
double delta4(int i, int j, int k, int l);
void Matrix_Inversion(double M[3][3], double Inv_M[3][3]);
double Determinant(double M[3][3]);
double box_volume(double a[3], double b[3], double c[3]);
double vec_dot(double a[3], double b[3]);
void vec_cross(double a[3], double b[3], double c[3]);

void make_Cijkl(double Cijkl[3][3][3][3], double Cij[6][6]);
void make_Sijkl(double Sijkl[3][3][3][3], double Sij[6][6]);

void make_Gik(double Cijkl[3][3][3][3], double r_vec[3], double Gir[3][3]);
void make_Gik_j(double Cijkl[3][3][3][3], double r_vec[3], double Gir_s[3][3][3]);
void make_Gik_jl(double Cijkl[3][3][3][3], double r_vec[3], double Gir_sm[3][3][3][3]);

double Eint_pair(double Cijkl[3][3][3][3], double r_vec[3], double Pij[3][3]);

#ifdef ORIENTATION
void make_Pij(double Cijkl[3][3][3][3], char **argv, double Omega1, double Omega2, double Pij[3][3]);
void read_input(double Cij[6][6], double Sij[6][6],
                double box_ref_1[3], double box_ref_2[3], double box_ref_3[3],
                double *a_lattice_ref, double *Omega1, double *Omega2);
#endif

#ifdef ABINITIO
void make_Pij(double Cijkl[3][3][3][3], double stress_ref[3][3], double stress_def[3][3],
              double box_ref_1[3], double box_ref_2[3], double box_ref_3[3],
              double box_def_1[3], double box_def_2[3], double box_def_3[3],
              double Pij[3][3]);

void read_input(double Cij[6][6], double Sij[6][6],
                double box_ref_1[3], double box_ref_2[3], double box_ref_3[3],
		double *a_lattice_ref, double stress_ref[3][3],
                double box_def_1[3], double box_def_2[3], double box_def_3[3],
                double *a_lattice_def, double stress_def[3][3]);
#endif

/********************************************************************************/

int main(int argc, char **argv){

    double Cij[6][6];
    double Sij[6][6];

    double box_ref_1[3];
    double box_ref_2[3];
    double box_ref_3[3];
    double a_lattice_ref;

    #ifdef ORIENTATION
    double Omega1;
    double Omega2;
    read_input(Cij, Sij, 
               box_ref_1, box_ref_2, box_ref_3, 
               &a_lattice_ref, &Omega1, &Omega2);
    #endif

    #ifdef ABINITIO
    double stress_ref[3][3];

    double box_def_1[3];
    double box_def_2[3];
    double box_def_3[3];
    double a_lattice_def;
    double stress_def[3][3];

    read_input(Cij, Sij, 
               box_ref_1, box_ref_2, box_ref_3, 
               &a_lattice_ref, stress_ref, 
               box_def_1, box_def_2, box_def_3, 
               &a_lattice_def, stress_def);
    #endif

    for (int i = 0; i < 3; ++i){
        box_ref_1[i] *= a_lattice_ref;
        box_ref_2[i] *= a_lattice_ref;
        box_ref_3[i] *= a_lattice_ref;
    }

    #ifdef ABINITIO
    for (int i = 0; i < 3; ++i){
        box_def_1[i] *= a_lattice_def;
        box_def_2[i] *= a_lattice_def;
        box_def_3[i] *= a_lattice_def;
    }
    #endif

    double Cijkl[3][3][3][3];
    double Sijkl[3][3][3][3];
    make_Cijkl(Cijkl, Cij); 
    make_Sijkl(Sijkl, Sij); 

    double P[3][3];  //dipole tensor

    #ifdef ORIENTATION
    make_Pij(Cijkl, argv, Omega1, Omega2, P);
    #endif
    #ifdef ABINITIO
    make_Pij(Cijkl, stress_ref, stress_def, box_ref_1, box_ref_2, box_ref_3, box_def_1, box_def_2, box_def_3, P);
    #endif

    cout << setiosflags(ios::scientific) << setprecision(16);
    cout << "Pij (eV) = " << "\n";
    cout << P[0][0] << " " << P[0][1] << " " << P[0][2] << "\n";
    cout << P[1][0] << " " << P[1][1] << " " << P[1][2] << "\n";
    cout << P[2][0] << " " << P[2][1] << " " << P[2][2] << "\n";
    cout << "\n";

    cout << "P11 = " << P[0][0] << " eV \n";
    cout << "P22 = " << P[1][1] << " eV \n";
    cout << "P33 = " << P[2][2] << " eV \n";
    cout << "P12 = " << P[0][1] << " eV \n";
    cout << "P23 = " << P[1][2] << " eV \n";
    cout << "P31 = " << P[2][0] << " eV \n";
    cout << "\n";

    #ifdef NOCORRECT
    int Range_Neigh = 30; // If defined NOCORRECT, one should not define SPERICALSUM. Using cubic sum can achieve 4 significant figures accuracy. 
    #else
    int Range_Neigh = 10; // Using 10 should be enough for spherical summation, but can be increased for checking.
    #endif

    int Range_x = Range_Neigh;
    int Range_y = Range_Neigh;
    int Range_z = Range_Neigh;
    int Range_r_sq = Range_Neigh*Range_Neigh;

    double Eint_DD = 0e0;
    for (int p = -Range_x; p <= Range_x; ++p){
        for (int q = -Range_y; q <= Range_y; ++q){
            for (int r = -Range_z; r <= Range_z; ++r){
                if (!(p == 0 && q == 0 && r == 0)){
                    #ifdef SPHERICALSUM
                    double r0_sq = p*p + q*q + r*r;
                    if (r0_sq < Range_r_sq+ 0.01){
                    #endif
                       double r_vec[3];
                       for (int m = 0 ; m < 3; ++m)
                           r_vec[m] = box_def_1[m]*p + box_def_2[m]*q + box_def_3[m]*r;
                       Eint_DD += Eint_pair(Cijkl, r_vec, P);
                    #ifdef SPHERICALSUM
                    }
                    #endif
                }
            }
        }
    }
    Eint_DD /= 2.0;

    double Eint_corr = 0e0;

    //Gaussian Quadrature: 9 points
    int Ngq = 9;
    double gq_pt[Ngq], gq_wt[Ngq];
    gq_pt[0] =  0.0000000000000000;
    gq_pt[1] =  0.8360311073266358;
    gq_pt[2] = -0.8360311073266358;
    gq_pt[3] =  0.9681602395076261;
    gq_pt[4] = -0.9681602395076261;
    gq_pt[5] =  0.3242534234038089;
    gq_pt[6] = -0.3242534234038089;
    gq_pt[7] =  0.6133714327005904;
    gq_pt[8] = -0.6133714327005904;
    gq_wt[0] =  0.3302393550012598;
    gq_wt[1] =  0.1806481606948574;
    gq_wt[2] =  0.1806481606948574;
    gq_wt[3] =  0.0812743883615744;
    gq_wt[4] =  0.0812743883615744;
    gq_wt[5] =  0.3123470770400029;
    gq_wt[6] =  0.3123470770400029;
    gq_wt[7] =  0.2606106964029354;
    gq_wt[8] =  0.2606106964029354;

    double Volume = box_volume(box_def_1, box_def_2, box_def_3);
    cout << "Defected box volume = " << Volume << " A^3\n";

    #ifndef NOCORRECT

    double int_Gik_jl[3][3][3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    int_Gik_jl[i][k][j][l] = 0e0;

    //We only tested the rectangular reference box. 
    //The correctness of the integration for orthorhombic cell is based on faith. 


    for (int p = -Range_x; p <= Range_x; ++p){
        for (int q = -Range_y; q <= Range_y; ++q){
            for (int r = -Range_z; r <= Range_z; ++r){
                if (!(p == 0 && q == 0 && r == 0)){
                    #ifdef SPHERICALSUM
                    double r0_sq = p*p + q*q + r*r;
                    if (r0_sq < Range_r_sq+ 0.01){
                    #endif
                       double r_vec[3];
                       for (int m = 0 ; m < 3; ++m)
                           r_vec[m] = box_def_1[m]*p + box_def_2[m]*q + box_def_3[m]*r;

                        for (int n_gq_1 = 0; n_gq_1 < Ngq; ++n_gq_1){
                            for (int n_gq_2 = 0; n_gq_2 < Ngq; ++n_gq_2){
                                for (int n_gq_3 = 0; n_gq_3 < Ngq; ++n_gq_3){
                                    double r_vec_0[3];
                                    r_vec_0[0] = r_vec[0] + (gq_pt[n_gq_1]*box_def_1[0] + gq_pt[n_gq_2]*box_def_2[0] + gq_pt[n_gq_3]*box_def_3[0])/2e0;
                                    r_vec_0[1] = r_vec[1] + (gq_pt[n_gq_1]*box_def_1[1] + gq_pt[n_gq_2]*box_def_2[1] + gq_pt[n_gq_3]*box_def_3[1])/2e0;
                                    r_vec_0[2] = r_vec[2] + (gq_pt[n_gq_1]*box_def_1[2] + gq_pt[n_gq_2]*box_def_2[2] + gq_pt[n_gq_3]*box_def_3[2])/2e0;
                                    double Gik_jl[3][3][3][3];
                                    make_Gik_jl(Cijkl, r_vec_0, Gik_jl);
                                    for (int i = 0; i < 3; ++i)
                                        for (int j = 0; j < 3; ++j)
                                            for (int k = 0; k < 3; ++k)
                                                for (int l = 0; l < 3; ++l)
                                                    int_Gik_jl[i][k][j][l] += gq_wt[n_gq_1]*gq_wt[n_gq_2]*gq_wt[n_gq_3]*Gik_jl[i][k][j][l];
                                }
                            }
                        }
                    #ifdef SPHERICALSUM
                    }
                    #endif
                }
            }
        }
    }

    double M[3][3];
    for (int i = 0; i < 3; ++i){
        M[0][i] = box_def_1[i];
        M[1][i] = box_def_2[i];
        M[2][i] = box_def_3[i];
    }   
    double Jacobian3D = Determinant(M)/8e0;
    // rx = (ux + 1)/2 * Lxx + (uy + 1)/2 * Lyx + (uz + 1)/2 * Lzx
    // ry = (ux + 1)/2 * Lxy + (uy + 1)/2 * Lyy + (uz + 1)/2 * Lzy
    // rz = (ux + 1)/2 * Lxz + (uy + 1)/2 * Lyz + (uz + 1)/2 * Lzz
    // Jacobian = | \parital (rx, ry, rz) / \partial (ux, uy, uz)|
    //          = Determinant(Box_volume_tensor) / 8
    cout << setprecision(16);
    cout << "Jacobian3D = " << Jacobian3D << " A^3\n\n";

    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            for (int j = 0; j < 3; ++j)
                for (int l = 0; l < 3; ++l)
                    int_Gik_jl[i][k][j][l] *= Jacobian3D/Volume;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    Eint_corr += P[i][j]*int_Gik_jl[i][k][j][l]*P[k][l];
    Eint_corr /= 2.0;

    #endif /*NOCORRECT*/


    double Estrain_corr = 0e0;

    // create 6 surface vectors
    double nn[6][3];
    vec_cross(box_def_1, box_def_2, nn[0]);    //top
    vec_cross(box_def_2, box_def_3, nn[1]);    //right
    vec_cross(box_def_3, box_def_1, nn[2]);    //back
    vec_cross(box_def_2, box_def_1, nn[3]);    //bottom
    vec_cross(box_def_3, box_def_2, nn[4]);    //left
    vec_cross(box_def_1, box_def_3, nn[5]);    //front
  
    double Jacobian2D[6] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    for (int n = 0; n < 6; ++n)
        Jacobian2D[n] = sqrt(vec_dot(nn[n],nn[n]))/4e0;


    //make them unit vector
    for (int n = 0; n < 6; ++n){
        double nn_norm = sqrt(nn[n][0]*nn[n][0] + nn[n][1]*nn[n][1] + nn[n][2]*nn[n][2]);
        for (int a = 0; a < 3; ++a){
            nn[n][a] /= nn_norm;
        }
    }

    double int_Gik_j_P_ka_n_a[6][3][3];  // 6 surface and i, j indexes
    for (int n = 0; n < 6; ++n)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                int_Gik_j_P_ka_n_a[n][i][j] = 0e0;

    for (int n = 0; n < 6; ++n){ //for 6 surfaces

        for (int n_gq_1 = 0; n_gq_1 < Ngq; ++n_gq_1){
            for (int n_gq_2 = 0; n_gq_2 < Ngq; ++n_gq_2){
                double r_vec_0[3] = {0e0, 0e0, 0e0};
                if (n == 0){  //top
                    r_vec_0[0] = (gq_pt[n_gq_1]*box_def_1[0] + gq_pt[n_gq_2]*box_def_2[0] + box_def_3[0])/2e0;
                    r_vec_0[1] = (gq_pt[n_gq_1]*box_def_1[1] + gq_pt[n_gq_2]*box_def_2[1] + box_def_3[1])/2e0;
                    r_vec_0[2] = (gq_pt[n_gq_1]*box_def_1[2] + gq_pt[n_gq_2]*box_def_2[2] + box_def_3[2])/2e0;
                } else if (n == 1) { //right
                    r_vec_0[0] = (box_def_1[0] + gq_pt[n_gq_1]*box_def_2[0] + gq_pt[n_gq_2]*box_def_3[0])/2e0;
                    r_vec_0[1] = (box_def_1[1] + gq_pt[n_gq_1]*box_def_2[1] + gq_pt[n_gq_2]*box_def_3[1])/2e0;
                    r_vec_0[2] = (box_def_1[2] + gq_pt[n_gq_1]*box_def_2[2] + gq_pt[n_gq_2]*box_def_3[2])/2e0;
                } else if (n == 2) { //back
                    r_vec_0[0] = (gq_pt[n_gq_1]*box_def_1[0] + box_def_2[0] + gq_pt[n_gq_2]*box_def_3[0])/2e0;
                    r_vec_0[1] = (gq_pt[n_gq_1]*box_def_1[1] + box_def_2[1] + gq_pt[n_gq_2]*box_def_3[1])/2e0;
                    r_vec_0[2] = (gq_pt[n_gq_1]*box_def_1[2] + box_def_2[2] + gq_pt[n_gq_2]*box_def_3[2])/2e0;
                } else if (n == 3){  //bottom
                    r_vec_0[0] = (gq_pt[n_gq_1]*box_def_1[0] + gq_pt[n_gq_2]*box_def_2[0] - box_def_3[0])/2e0;
                    r_vec_0[1] = (gq_pt[n_gq_1]*box_def_1[1] + gq_pt[n_gq_2]*box_def_2[1] - box_def_3[1])/2e0;
                    r_vec_0[2] = (gq_pt[n_gq_1]*box_def_1[2] + gq_pt[n_gq_2]*box_def_2[2] - box_def_3[2])/2e0;
                } else if (n == 4) { //left
                    r_vec_0[0] = (-box_def_1[0] + gq_pt[n_gq_1]*box_def_2[0] + gq_pt[n_gq_2]*box_def_3[0])/2e0;
                    r_vec_0[1] = (-box_def_1[1] + gq_pt[n_gq_1]*box_def_2[1] + gq_pt[n_gq_2]*box_def_3[1])/2e0;
                    r_vec_0[2] = (-box_def_1[2] + gq_pt[n_gq_1]*box_def_2[2] + gq_pt[n_gq_2]*box_def_3[2])/2e0;
                } else if (n == 5) { //front
                    r_vec_0[0] = (gq_pt[n_gq_1]*box_def_1[0] - box_def_2[0] + gq_pt[n_gq_2]*box_def_3[0])/2e0;
                    r_vec_0[1] = (gq_pt[n_gq_1]*box_def_1[1] - box_def_2[1] + gq_pt[n_gq_2]*box_def_3[1])/2e0;
                    r_vec_0[2] = (gq_pt[n_gq_1]*box_def_1[2] - box_def_2[2] + gq_pt[n_gq_2]*box_def_3[2])/2e0;
                }
                double Gik_j[3][3][3];
                make_Gik_j(Cijkl, r_vec_0, Gik_j);

                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        for (int k = 0; k < 3; ++k)
                            for (int a = 0; a < 3; ++a)
                                int_Gik_j_P_ka_n_a[n][i][j] += gq_wt[n_gq_1]*gq_wt[n_gq_2]*Gik_j[i][k][j]*P[k][a]*nn[n][a];

            }
        }
    }




    double epsilon_corr[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            epsilon_corr[i][j] = 0e0;

    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
             for (int n = 0; n < 6; ++n){
                int_Gik_j_P_ka_n_a[n][i][j] *= Jacobian2D[n]; //for gaussian quadrature
                epsilon_corr[i][j] += int_Gik_j_P_ka_n_a[n][i][j];
            }
            epsilon_corr[i][j] /= Volume;
            //cout << epsilon_corr[i][j] << " " << P[i][j] << '\n';
            Estrain_corr += epsilon_corr[i][j]*P[i][j];
        }
   }

    Estrain_corr /= 2.0; 

    double Eint = Eint_DD - Eint_corr;
    double Eel_corr = Eint - Estrain_corr;

    cout << "Eint_DD       = " <<  Eint_DD	  <<  " eV\n";
    cout << "Eint_corr     = " <<  Eint_corr      <<  " eV\n";
    cout << "Estrain_corr  = " <<  Estrain_corr      <<  " eV\n";
    cout << "Eint          = Eint_DD - Eint_corr   = " << Eint << " eV\n";
    cout << "Eel_corr      = Eint - Estrain_corr   = " << Eel_corr << " eV\n";
    cout << "\n";
 
    cout << "Eint_DD      = 1/2 * Sum Eint_pair (R_n)\n";
    cout << "Eint_corr    = 1/2 int_V_cell Sum Eint_pair (R_n - r) dV\n";
    cout << "Estrain_corr = 1/2 int_V_cell Pij int_S Gik_j (R) P_ka n_a dS\n";
    cout << "\n";

    cout << "E_F(isolated defect) = E_F(defect with periodic images) - Eel_corr\n";
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

    cout << "Relaxation Volume Tensor Omega_ij (A^3) = \n";
    cout << Omega[0][0] << " " << Omega[0][1] << " " << Omega[0][2] << "\n";
    cout << Omega[1][0] << " " << Omega[1][1] << " " << Omega[1][2] << "\n";
    cout << Omega[2][0] << " " << Omega[2][1] << " " << Omega[2][2] << "\n";
    cout << "\n";

    cout << "Relaxation Volume = " << Omega[0][0] + Omega[1][1] + Omega[2][2] << " A^3 \n";
    cout << "\n";

    cout << "Omega11 = " << Omega[0][0] << " A^3 \n";
    cout << "Omega22 = " << Omega[1][1] << " A^3 \n";
    cout << "Omega33 = " << Omega[2][2] << " A^3 \n";
    cout << "Omega12 = " << Omega[0][1] << " A^3 \n";
    cout << "Omega23 = " << Omega[1][2] << " A^3 \n";
    cout << "Omega31 = " << Omega[2][0] << " A^3 \n";
    cout << "\n";

}            
               
double delta(int i, int j){

    if (i == j){
        return 1.0;
    } else {
        return 0.0;
    }
}

double  delta4(int i, int j, int k, int l){

    if (i == j && i == k &&  i == l){
        return 1.0;
    } else {
        return 0.0;
    }
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

void make_Gik(double Cijkl[3][3][3][3], double r_vec[3], double Gir[3][3]){

    //Gik is calculated according to D. M. Barnett Phys. Stat. Sol. (b) 49, 741 (1972)

    double Pi = 3.141592653589793;

    double r0 = sqrt(pow(r_vec[0],2) + pow(r_vec[1],2) + pow(r_vec[2],2));

    double T[3];
    for (int i = 0; i < 3; ++i) T[i] = r_vec[i]/r0;

    //Note: We followed the spherical coordinate system that is used in Barnett's paper. 

    double cosPhi, sinPhi, cosTheta, sinTheta;
    if (fabs(T[0]) < 1e-8 && fabs(T[1]) < 1e-8){
        cosPhi = 1e0;
        sinPhi = 0e0;
        cosTheta = 1e0;
        sinTheta = 0e0;
    } else {
        cosPhi = T[2];
        sinPhi = sin(acos(T[2]));
        cosTheta=T[0]/sinPhi;
        sinTheta=T[1]/sinPhi;
    }

    double a[3], b[3];
    a[0] = sinTheta;
    a[1] = -cosTheta;
    a[2] = 0e0;
    b[0] = cosPhi*cosTheta;
    b[1] = cosPhi*sinTheta;
    b[2] = -sinPhi;

    int NPsi = 20;
    //Note: Using NPsi = 20 should give answer up to 4 significant figures correct.
    double dPsi = Pi/NPsi;

    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
            Gir[i][r] = 0e0;

    for (int nn = 0; nn < NPsi; ++nn){ //perform integration

        double Psi = Pi/NPsi*(nn+0.5);
        double cosPsi = cos(Psi);
        double sinPsi = sin(Psi);

        double z[3];
        for (int i = 0; i < 3; ++i) z[i] = a[i]*cosPsi + b[i]*sinPsi;

        double M[3][3], Inv_M[3][3];

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                M[i][j] = 0e0;

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int r = 0; r < 3; ++r)
                    for (int s = 0; s < 3; ++s)
                        if (Cijkl[i][j][r][s] > 0e0) M[i][r] += Cijkl[i][j][r][s]*z[j]*z[s];

        Matrix_Inversion(M, Inv_M);

        for (int i = 0; i < 3; ++i)
            for (int r = 0; r < 3; ++r)
                Gir[i][r] += Inv_M[i][r];
    }

    double factor = dPsi/(4e0*pow(Pi,2)*r0);
    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
            Gir[i][r] *= factor;
    
}

void make_Gik_j(double Cijkl[3][3][3][3], double r_vec[3], double Gir_s[3][3][3]){

    //Gik_j is calculated according to D. M. Barnett Phys. Stat. Sol. (b) 49, 741 (1972)

    double Pi = 3.141592653589793;

    double rsq = pow(r_vec[0],2) + pow(r_vec[1],2) + pow(r_vec[2],2);
    double r0 = sqrt(rsq);

    double T[3];
    for (int i = 0; i < 3; ++i) T[i] = r_vec[i]/r0;

    //Note: We followed the spherical coordinate system that is used in Barnett's paper. 

    double cosPhi, sinPhi, cosTheta, sinTheta;
    if (fabs(T[0]) < 1e-8 && fabs(T[1]) < 1e-8){
        cosPhi = 1e0;
        sinPhi = 0e0;
        cosTheta = 1e0;
        sinTheta = 0e0;
    } else {
        cosPhi = T[2];
        sinPhi = sin(acos(T[2]));
        cosTheta=T[0]/sinPhi;
        sinTheta=T[1]/sinPhi;
    }

    double a[3], b[3];
    a[0] = sinTheta;
    a[1] = -cosTheta;
    a[2] = 0e0;
    b[0] = cosPhi*cosTheta;
    b[1] = cosPhi*sinTheta;
    b[2] = -sinPhi;

    int NPsi = 20;
    //Note: Using NPsi = 20 should give answer up to 4 significant figures correct.
    double dPsi = Pi/NPsi;

    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
            for (int s = 0; s < 3; ++s)
                    Gir_s[i][r][s] = 0e0;

    for (int nn = 0; nn < NPsi; ++nn){ //perform integration

        double Psi = Pi/NPsi*(nn+0.5);
        double cosPsi = cos(Psi);
        double sinPsi = sin(Psi);

        double z[3];
        for (int i = 0; i < 3; ++i) z[i] = a[i]*cosPsi + b[i]*sinPsi;

        double M[3][3], Inv_M[3][3], F[3][3];

        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                M[i][j] = 0e0;
                F[i][j] = 0e0;
            }
        }

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int r = 0; r < 3; ++r)
                    for (int s = 0; s < 3; ++s)
                        if (Cijkl[i][j][r][s] > 0e0) M[i][r] += Cijkl[i][j][r][s]*z[j]*z[s];

        Matrix_Inversion(M, Inv_M);

        for (int i = 0; i < 3; ++i)
            for (int r = 0; r < 3; ++r)
                for (int j = 0; j < 3; ++j)
                    for (int p = 0; p < 3; ++p)
                        for (int n = 0; n < 3; ++n)
                            for (int w = 0; w < 3; ++w)
                                if (Cijkl[j][p][n][w] > 0e0) 
                                    F[i][r] += Cijkl[j][p][n][w]*Inv_M[i][j]*Inv_M[n][r]*(z[p]*T[w] + z[w]*T[p]);

        for (int i = 0; i < 3; ++i)
            for (int r = 0; r < 3; ++r)
                for (int s = 0; s < 3; ++s)
                        Gir_s[i][r][s] += -T[s]*Inv_M[i][r] + z[s]*F[i][r];
    }

    double factor = dPsi/(4e0*pow(Pi,2)*rsq);
    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
            for (int s = 0; s < 3; ++s)
                    Gir_s[i][r][s] *= factor;
    
}

void make_Gik_jl(double Cijkl[3][3][3][3], double r_vec[3], double Gir_sm[3][3][3][3]){

    //Gik_jl is calculated according to D. M. Barnett Phys. Stat. Sol. (b) 49, 741 (1972)

    double Pi = 3.141592653589793;

    double r0 = sqrt(pow(r_vec[0],2) + pow(r_vec[1],2) + pow(r_vec[2],2));

    double T[3];
    for (int i = 0; i < 3; ++i) T[i] = r_vec[i]/r0;

    //Note: We followed the spherical coordinate system that is used in Barnett's paper. 

    double cosPhi, sinPhi, cosTheta, sinTheta;
    if (fabs(T[0]) < 1e-8 && fabs(T[1]) < 1e-8){
        cosPhi = 1e0;
        sinPhi = 0e0;
        cosTheta = 1e0;
        sinTheta = 0e0;
    } else {
        cosPhi = T[2];
        sinPhi = sin(acos(T[2]));
        cosTheta=T[0]/sinPhi;
        sinTheta=T[1]/sinPhi;
    }

    double a[3], b[3];
    a[0] = sinTheta;
    a[1] = -cosTheta;
    a[2] = 0e0;
    b[0] = cosPhi*cosTheta;
    b[1] = cosPhi*sinTheta;
    b[2] = -sinPhi;

    int NPsi = 20; 
    //Note: Using NPsi = 20 should give answer up to 4 significant figures correct.
    double dPsi = Pi/NPsi;

    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
            for (int s = 0; s < 3; ++s)
                for (int m = 0; m < 3; ++m)
                    Gir_sm[i][r][s][m] = 0e0;

    for (int nn = 0; nn < NPsi; ++nn){ //perform integration

        double Psi = Pi/NPsi*(nn+0.5);
        double cosPsi = cos(Psi);
        double sinPsi = sin(Psi);

        double z[3];
        for (int i = 0; i < 3; ++i) z[i] = a[i]*cosPsi + b[i]*sinPsi;

        double M[3][3], Inv_M[3][3], A[3][3], F[3][3];

        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                M[i][j] = 0e0;
                A[i][j] = 0e0;
                F[i][j] = 0e0;
            }
        }

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int r = 0; r < 3; ++r)
                    for (int s = 0; s < 3; ++s)
                        if (Cijkl[i][j][r][s] > 0e0) M[i][r] += Cijkl[i][j][r][s]*z[j]*z[s];

        Matrix_Inversion(M, Inv_M);

        for (int i = 0; i < 3; ++i)
            for (int r = 0; r < 3; ++r)
                for (int j = 0; j < 3; ++j)
                    for (int p = 0; p < 3; ++p)
                        for (int n = 0; n < 3; ++n)
                            for (int w = 0; w < 3; ++w)
                                if (Cijkl[j][p][n][w] > 0e0) 
                                    F[i][r] += Cijkl[j][p][n][w]*Inv_M[i][j]*Inv_M[n][r]*(z[p]*T[w] + z[w]*T[p]);


        for (int i = 0; i < 3; ++i)
            for (int r = 0; r < 3; ++r)
                for (int j = 0; j < 3; ++j)
                    for (int p = 0; p < 3; ++p)
                        for (int n = 0; n < 3; ++n)
                            for (int w = 0; w < 3; ++w)
                                if (Cijkl[j][p][n][w] > 0e0) 
                                    A[i][r] += Cijkl[j][p][n][w]*((z[p]*T[w]+z[w]*T[p])*(F[i][j]*Inv_M[n][r]+Inv_M[i][j]*F[n][r])
                                               -2e0*Inv_M[i][j]*Inv_M[n][r]*T[p]*T[w]);

        for (int i = 0; i < 3; ++i)
            for (int r = 0; r < 3; ++r)
                for (int s = 0; s < 3; ++s)
                    for (int m = 0; m < 3; ++m)
                        Gir_sm[i][r][s][m] += 2e0*T[s]*T[m]*Inv_M[i][r] - 2e0*(z[s]*T[m]+z[m]*T[s])*F[i][r] + z[s]*z[m]*A[i][r];
    }

    double factor = dPsi/(4e0*pow(Pi,2)*pow(r0,3));
    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
            for (int s = 0; s < 3; ++s)
                for (int m = 0; m < 3; ++m)
                    Gir_sm[i][r][s][m] *= factor;
    
}

void Matrix_Inversion(double M[3][3], double Inv_M[3][3]){

    Inv_M[0][0] = M[1][1]*M[2][2] - M[2][1]*M[1][2];
    Inv_M[0][1] = M[0][2]*M[2][1] - M[0][1]*M[2][2];
    Inv_M[0][2] = M[0][1]*M[1][2] - M[0][2]*M[1][1];
    Inv_M[1][0] = M[1][2]*M[2][0] - M[1][0]*M[2][2];
    Inv_M[1][1] = M[0][0]*M[2][2] - M[0][2]*M[2][0];
    Inv_M[1][2] = M[1][0]*M[0][2] - M[0][0]*M[1][2];
    Inv_M[2][0] = M[1][0]*M[2][1] - M[2][0]*M[1][1];
    Inv_M[2][1] = M[2][0]*M[0][1] - M[0][0]*M[2][1];
    Inv_M[2][2] = M[0][0]*M[1][1] - M[1][0]*M[0][1];

    double Det = M[0][0]*Inv_M[0][0]
               + M[0][1]*Inv_M[1][0]
               + M[0][2]*Inv_M[2][0];

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            Inv_M[i][j] /= Det;
}

double Determinant(double M[3][3]){

    double Inv_M[3][3];
    Inv_M[0][0] = M[1][1]*M[2][2] - M[2][1]*M[1][2];
    Inv_M[1][0] = M[1][2]*M[2][0] - M[1][0]*M[2][2];
    Inv_M[2][0] = M[1][0]*M[2][1] - M[2][0]*M[1][1];

    double Det = M[0][0]*Inv_M[0][0]
               + M[0][1]*Inv_M[1][0]
               + M[0][2]*Inv_M[2][0];

    return Det;
}

double box_volume(double a[3], double b[3], double c[3]){

    double d[3];
    d[0] = a[1]*b[2] - b[1]*a[2];
    d[1] = a[2]*b[0] - b[2]*a[0];
    d[2] = a[0]*b[1] - b[0]*a[1];

    double volume;
    volume = c[0]*d[0] + c[1]*d[1] + c[2]*d[2];
    volume = fabs(volume);

    return volume;

}

double vec_dot(double a[3], double b[3]){

    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void vec_cross(double a[3], double b[3], double c[3]){

    c[0] = a[1]*b[2] - b[1]*a[2];
    c[1] = a[2]*b[0] - b[2]*a[0];
    c[2] = a[0]*b[1] - b[0]*a[1];

}


double Eint_pair(double Cijkl[3][3][3][3], double r_vec[3], double Pij[3][3]){

    double G[3][3][3][3];
    make_Gik_jl(Cijkl, r_vec, G);

    double E_int = 0e0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l){
                    E_int += Pij[i][j]*Pij[k][l]*G[i][k][j][l];
                }
    return E_int;
}

#ifdef ORIENTATION
void make_Pij(double Cijkl[3][3][3][3], char **argv, double Omega1, double Omega2, double Pij[3][3]){

    const double Pi = 3.141592653589793;

    double theta0 = atof(argv[1]);
    double phi0   = atof(argv[2]);

    double theta = theta0/180e0*Pi;
    double phi   = phi0/180e0*Pi;

    double unit_n[3]; // the orientation of defect
    unit_n[0] = cos(phi)*sin(theta);
    unit_n[1] = sin(phi)*sin(theta);
    unit_n[2] = cos(theta);

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            Pij[i][j] = 0e0;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    if (Cijkl[i][j][k][l] > 0e0) 
                        Pij[i][j] += Cijkl[i][j][k][l]*(Omega1*unit_n[k]*unit_n[l] + Omega2/3e0*delta(k,l));

}
#endif

#ifdef ABINITIO
void make_Pij(double Cijkl[3][3][3][3], double stress_ref[3][3], double stress_def[3][3],
              double box_ref_1[3], double box_ref_2[3], double box_ref_3[3],
              double box_def_1[3], double box_def_2[3], double box_def_3[3],
              double Pij[3][3]){

    double delta_box[3][3], box_ref[3][3], Inv_box_ref[3][3];

    for (int i = 0; i < 3; ++i){
        delta_box[0][i] = box_def_1[i] - box_ref_1[i];
        delta_box[1][i] = box_def_2[i] - box_ref_2[i];
        delta_box[2][i] = box_def_3[i] - box_ref_3[i];
        box_ref[0][i] = box_ref_1[i];
        box_ref[1][i] = box_ref_2[i];
        box_ref[2][i] = box_ref_3[i];
    }

    Matrix_Inversion(box_ref, Inv_box_ref);
   
    double strain_apply[3][3], stress_apply[3][3];  

    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
             strain_apply[i][j] = 0e0;
             stress_apply[i][j] = 0e0;
        }
    }
 
    //calculate the applied strain from the deformation of defect simulation box   
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                strain_apply[i][j] += (Inv_box_ref[i][k]*delta_box[k][j] 
                                     + delta_box[i][k]*Inv_box_ref[k][j])/2e0;
           

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    stress_apply[i][j] += Cijkl[i][j][k][l]*strain_apply[k][l];

    double volume = box_volume(box_ref_1, box_ref_2, box_ref_3);
    
    double eVperAcubic_to_GPa = 160.21766208;
    double stress[3][3];

    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            stress[i][j] = stress_ref[i][j] - stress_def[i][j];
            #ifdef STRESSeV
            //If stresses is total stress in eV,
            //there is no need to be multiplied by the volume,
            //because it was multiplied in the DFT program (VASP).
            Pij[i][j] = volume*stress_apply[i][j] - stress[i][j];
            #endif
            #ifdef STRESSGPa
            //Stresses in GPa
            //Need to be converted to eV A^(-3) 
            //and be multiplied by volume in A^3.
            stress[i][j] /= eVperAcubic_to_GPa;         
            Pij[i][j] = volume*(stress_apply[i][j] - stress[i][j]);
            #endif
        } 
    }

    double E_strain_apply = 0e0;
    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            E_strain_apply += strain_apply[i][j]*(0.5*volume*stress_apply[i][j] - Pij[i][j]);
        }
    }

    cout << setprecision(16);
    cout << "Energy due to the deformation of defect simulation box\n";
    cout << "E_strain_apply = 0.5*V_ref*(epsilon_ij C_ijkl epsilon_kl) - Pij epsilon_ij = " << E_strain_apply << " eV\n";
    cout << "Reference box volume = " << volume << " A^3\n";
    cout << "\n";
}
#endif

#ifdef ORIENTATION
void read_input(double Cij[6][6], double Sij[6][6],
                double box_ref_1[3], double box_ref_2[3], double box_ref_3[3],
                double *a_lattice_ref, double *Omega1, double *Omega2)                      
#endif
#ifdef ABINITIO
void read_input(double Cij[6][6], double Sij[6][6],
                double box_ref_1[3], double box_ref_2[3], double box_ref_3[3],
                double *a_lattice_ref, double stress_ref[3][3], 
                double box_def_1[3], double box_def_2[3], double box_def_3[3],
                double *a_lattice_def, double stress_def[3][3])
#endif
{
    ifstream in_file_1("input_data");
    ifstream in_file_2("input_elastic");


    if (in_file_1) {
        cout << "Start reading input_data ... \n\n";
    } else {
        cout << "ERROR: You need to have an input_data file \n\n";
        exit(1);
    }


    while (!in_file_1.eof()){

        char variable[256];
        memset(variable, 0, 256);

        char value[256];
        memset(value, 0, 256);

        in_file_1 >> variable >> value;
        if (strcmp(variable, "box_ref_11") == 0) box_ref_1[0] = atof(value);
        if (strcmp(variable, "box_ref_12") == 0) box_ref_1[1] = atof(value);
        if (strcmp(variable, "box_ref_13") == 0) box_ref_1[2] = atof(value);
        if (strcmp(variable, "box_ref_21") == 0) box_ref_2[0] = atof(value);
        if (strcmp(variable, "box_ref_22") == 0) box_ref_2[1] = atof(value);
        if (strcmp(variable, "box_ref_23") == 0) box_ref_2[2] = atof(value);
        if (strcmp(variable, "box_ref_31") == 0) box_ref_3[0] = atof(value);
        if (strcmp(variable, "box_ref_32") == 0) box_ref_3[1] = atof(value);
        if (strcmp(variable, "box_ref_33") == 0) box_ref_3[2] = atof(value);
        if (strcmp(variable, "a_lattice_ref") == 0) *a_lattice_ref = atof(value);

        #ifdef ORIENTATION
        if (strcmp(variable, "Omega1") == 0) *Omega1 = atof(value);
        if (strcmp(variable, "Omega2") == 0) *Omega2 = atof(value);
        #endif

        #ifdef ABINITIO
        if (strcmp(variable, "box_def_11") == 0) box_def_1[0] = atof(value);
        if (strcmp(variable, "box_def_12") == 0) box_def_1[1] = atof(value);
        if (strcmp(variable, "box_def_13") == 0) box_def_1[2] = atof(value);
        if (strcmp(variable, "box_def_21") == 0) box_def_2[0] = atof(value);
        if (strcmp(variable, "box_def_22") == 0) box_def_2[1] = atof(value);
        if (strcmp(variable, "box_def_23") == 0) box_def_2[2] = atof(value);
        if (strcmp(variable, "box_def_31") == 0) box_def_3[0] = atof(value);
        if (strcmp(variable, "box_def_32") == 0) box_def_3[1] = atof(value);
        if (strcmp(variable, "box_def_33") == 0) box_def_3[2] = atof(value);
        if (strcmp(variable, "a_lattice_def") == 0) *a_lattice_def = atof(value);

        if (strcmp(variable, "stress11_ref") == 0) stress_ref[0][0] = atof(value);
        if (strcmp(variable, "stress12_ref") == 0) stress_ref[0][1] = atof(value);
        if (strcmp(variable, "stress13_ref") == 0) stress_ref[0][2] = atof(value);
        if (strcmp(variable, "stress21_ref") == 0) stress_ref[1][0] = atof(value);
        if (strcmp(variable, "stress22_ref") == 0) stress_ref[1][1] = atof(value);
        if (strcmp(variable, "stress23_ref") == 0) stress_ref[1][2] = atof(value);
        if (strcmp(variable, "stress31_ref") == 0) stress_ref[2][0] = atof(value);
        if (strcmp(variable, "stress32_ref") == 0) stress_ref[2][1] = atof(value);
        if (strcmp(variable, "stress33_ref") == 0) stress_ref[2][2] = atof(value);

        if (strcmp(variable, "stress11_def") == 0) stress_def[0][0] = atof(value);
        if (strcmp(variable, "stress12_def") == 0) stress_def[0][1] = atof(value);
        if (strcmp(variable, "stress13_def") == 0) stress_def[0][2] = atof(value);
        if (strcmp(variable, "stress21_def") == 0) stress_def[1][0] = atof(value);
        if (strcmp(variable, "stress22_def") == 0) stress_def[1][1] = atof(value);
        if (strcmp(variable, "stress23_def") == 0) stress_def[1][2] = atof(value);
        if (strcmp(variable, "stress31_def") == 0) stress_def[2][0] = atof(value);
        if (strcmp(variable, "stress32_def") == 0) stress_def[2][1] = atof(value);
        if (strcmp(variable, "stress33_def") == 0) stress_def[2][2] = atof(value);
        #endif
    }
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

