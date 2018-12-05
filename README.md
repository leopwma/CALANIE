# CALANIE

   Copyright (C) 2018 Culham Centre for Fusion Energy,
   United Kingdom Atomic Energy Authority, Oxfordshire OX14 3DB, UK

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

*******************************************************************************
   Program: CALANIE
   
            CALculation of ANIsotropic Elastic interaction energy of
            a defect in periodic boundary conditions
   
   Version: 1.9
   
   Date:    05 Dec 2018
   
   Author:  Pui-Wai (Leo) MA
   
   Contact: leo.ma@ukaea.uk
   
   Address: Culham Centre for Fusion Energy, OX14 3DB, United Kingdom

********************************************************************************
 
  This program calculate the elastic interaction energy of a defect with its
  images in periodc boundary conditions. It also output the elastic dipole 
  tensor and relaxation volume tensor.
  
  It requires two input files: input_data & input_elastic.
  "input_data" will be discussed here.

  "input_elastic" contains the information of elastic constant tensor.
  Format of it is mentioned in "make_input_elastic.py".

# Compilation:

  One can compile the program in 2 ways.

  1) g++ -DORIENTATION CALANIE_version.cpp

     This requires input of 
     
     (1) Reference simulation box vectors
     
     (2) Lattice unit for reference box vectors
     
     (3) Relaxation volumes Omega1 and Omega2
     
     (4) Orientation of the defect

     Put (1) to (3) in a file with name "input_data", and also the 
     "input_elastic". 

     Run the program by
     
     ./a.out theta phi
     
     where theta and phi are (4) Orientation of the defect

     Detail is in S. L. Dudarev and Pui-Wai Ma, 
     Physical Review Materials 2 (3), 033602 (2018)

  Or

  2) g++ -DABINITIO -DSTRESSeV CALANIE_version.cpp
     
     or
     
     g++ -DABINITIO -DSTRESSGPa CALANIE_version.cpp

     This require input of
     
     (1) Reference simulation box vectors
     
     (2) Lattice units for reference simulation box vectors
     
     (3) Total marcoscopic stress of the reference simulation box (in eV or GPa) 
     
     (4) Defect simulation box vectors
     
     (5) Lattice units for defect simulation box vectors
     
     (6) Total marcoscopic stress of the defect simulation box (in eV or GPa)

     Put (1) to (6) in a file with name "input_data"

     Run the program by
     
     ./a.out

     Note:
     
     -DSTRESSeV means total stresses multiply volume. Unit is in eV. It is 
     exactly the line of "Total" in vasp "OUTCAR". 

     -DSTRESSGPa means total stresses. Unit is in GPa. It is corresponding to 
     the line of "kB" in vasp "OUTCAR", where "kB" means kbar = 0.1GPa.

*******************************************************************************
 Optional:

(1) In the correct term, we are doing cubic sum. We can change it into
        spherical sum by adding -DSPHERICALSUM . 

(2) We observed that if we use cubic sum, the correction term tends to
        zero, if we consider more neighbour cells. One can ignore the 
        calculation of the correction term by adding -DNOCORRECT. The default 
        value of Range_Neigh changes from 10 to 30. However, it is much faster.  

(3) <b>(Very important!!!)</b> It is strongly recommended to add either 
        -DSPHERICALSUM or -DNOCORRECT to the compilation, but not both.
        
        For accuracy, use -DSPHERICALSUM.
        
        For speed, use -DNOCORRECT.

*******************************************************************************

  A sample "input_data" is given for both cases. "input_data_1" is for 1st 
  case. "input_data_2" is for 2nd case using -DSTRESSeV. 

  Sign convention of total marcoscopic stress follow VASP. 
  i.e positive means the simulation cell intends to expand.

