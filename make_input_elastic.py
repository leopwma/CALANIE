################################################################################
#
#   Copyright (C) 2018 Culham Centre for Fusion Energy,
#   United Kingdom Atomic Energy Authority, Oxfordshire OX14 3DB, UK
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
################################################################################
#
#   Program: make_input_elastic.py
#   Version: 1.0
#   Date:    June 2018
#   Author:  Pui-Wai (Leo) MA
#   Contact: leo.ma@ukaea.uk
#   Address: Culham Centre for Fusion Energy, OX14 3DB, United Kingdom
#
################################################################################
# 
#   This is a small program making the compliance matrix Sij = Cij^-1.
#   It produces one of the two input files for CALANIE, i.e. input_elastic.
#   It works with python2 and python3, with numpy. 
#   
#   It is not necessary to use this python program to produce "input_elastic",
#   if one can produce it according to the format mentioned below.
#
################################################################################
#  
#   The input file of this python program is "input_elastic_Cij".
#   Format:
#
#   line 1 --> comment
#   line 2 --> comment
#   line 6 to 11 --> C11 C12 C13 C14 C15 C16
#                    C21 C22 C23 C24 C25 C26
#                    C31 C32 C33 C34 C35 C36
#                    C41 C42 C43 C44 C45 C46
#                    C51 C52 C53 C54 C55 C56
#                    C61 C62 C63 C64 C65 C66
#   Cij is in unit of GPa.

#   The output file of this python program is "input_elastic".
#   
#   The first 1 to 11 lines are the same as "input_elastic_Cij"
#   line 12 to 17 --> S11 S12 S13 S14 S15 S16
#                     S21 S22 S23 S24 S25 S26
#                     S31 S32 S33 S34 S35 S36
#                     S41 S42 S43 S44 S45 S46
#                     S51 S52 S53 S54 S55 S56
#                     S61 S62 S63 S64 S65 S66
#   Sij is in unit of GPa^-1.
#
################################################################################

import numpy as np
from numpy.linalg import inv

Cij = np.zeros((6,6))
Sij = np.zeros((6,6))

f1 = open("input_elastic_Cij","r")
f2 = open("input_elastic","w")
for i, line in enumerate(f1):
    if i <=1:
        f2.write(line)
    if i > 1:
        x = line.split()
        x = np.array(x,dtype=float)
        Cij[i-2] = x

f1.close()

Sij = inv(Cij)

for i in range(6):
    f2.write("%23.16e %23.16e %23.16e %23.16e %23.16e %23.16e \n" % (Cij[i,0], Cij[i,1], Cij[i,2], Cij[i,3], Cij[i,4], Cij[i,5]))

for i in range(6):
    f2.write("%23.16e %23.16e %23.16e %23.16e %23.16e %23.16e \n" % (Sij[i,0], Sij[i,1], Sij[i,2], Sij[i,3], Sij[i,4], Sij[i,5]))

f2.close()




