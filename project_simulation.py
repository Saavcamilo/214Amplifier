#!/usr/bin/python

import random
import collections
from math import *
import sys
from math import *


Vdd = 2.5
Vss = -2.5
Imin = 1000
lam = .1
gam = 0.6
phi = 0.4
V_th = 0.5
uCox = 50e-6
uCox_p = 25e-6
outputCurrent, outputVoltages, outputParams, outputPerform = [], [], [], []

Cin = 100 * 10 ** -15
# Capacitance constants
Cgate = 2.03 * 10 ** -15
# Ratio of parasitic capacitance
Cdb = .33
# Assume Csb = Cdb 
Csb = .33
Cgd = .25

# Sweep pvbias range:
Vbias_min = -1.9
Vbias_max = -1.4
Vbias_inc = .00001
Vbias = Vbias_min

Vbiasp_min = 1.4
Vbiasp_max = 1.9
Vbiasp = Vbiasp_min#Vbiasp_min

RU = 40000.
RD = 60000.

W_L_1L = 40.
W_L_1bias = 32.
W_L_1 = 13.

W_L_2L = 20.
W_L_2bias = 32.
W_L_2 = 13.

W_L_3 = 13.
W_L_3bias = 32.
Amax = -1000.
# Loop through Vbiasp, determine current down
while(Vbiasp < Vbiasp_max):
    vovLoad = Vdd - Vbiasp - .5
    IdL = .5 * uCox_p * W_L_1L * vovLoad ** 2
    # RoL is given by IdL
    roL = 1/(lam * IdL)
    #Loop through bottom biasing
    while(Vbias < Vbias_max):
        vovB = Vbias - Vss - V_th    
        IdB1 = .5 * uCox * W_L_1bias * vovB ** 2
        vov1 = sqrt(2 * IdB1 / (W_L_1 * uCox))       
        roB = 1/(lam * IdL)
        # Body Effect Time
        # This equation is a taylor expansion around vov1 = 0
        Vsb = 1.60598 - .837037 * vov1
        Vs1 = Vsb + Vss
        Ires = IdL - IdB1
        # Vd1 = Current through it * R tot
        Vd1 = (Ires + Vdd/RU + Vss/RD)/(1/RU + 1/RD)
        if Vd1 > 2.5 or Vd1 < 0:
            Vbias += Vbias_inc
            continue
        vov1 = -Vs1-.5
        gmL = (2*IdL)/vovLoad
        gm1 = (2*IdB1)/(vov1)
        # Use vovPmos and Itop to determine W of M_L1
        Iu = (Vdd - Vd1)/RU
        Id = (Vd1 - Vss)/RD
        # Rtot = roL//Ru//Rd
        Rtot = roL * RU * RD / (roL * RU + roL * RD + RU * RD)
        A1 = -Rtot
        # Bandwidth Calculations for Stage 1
        # Calculate Cgs1 from tech. parameters
        Cgs1 = W_L_1 * Cgate
        # Approximate other caps as ratio with Cgs1 for simplicity
        Cdb1 = Cgs1 * Cdb
        Csb1 = Cgs1 * Csb
        Cgd1 = Cgs1 * Cgd
        #rMiller to separate common gate ro and place in parallel
        # Miller Transform -> input = Z/(1-K), output = Z * 1/(1-1/K)
        rMiller = roB / (1 - Vd1/ Vs1)
        R3db = roB/(1 + gm1*roB)
        Cgs1B = W_L_1bias * Cgate
        Cdb1B = Cgs1B * Cdb
        Cgd1B = Cgs1B * Cgd
        CtotB = Cgd1B + Cdb1B 
        Ctot = Cgs1 + Csb1 + CtotB + Cin
        #print 1/(Ctot * R3db * 2 * pi * 10 ** 6)
        # Approximate other caps as ratio with Cgs1 for simplicity

        Cgs1L = W_L_1L * Cgate
        Cdb1L = Cgs1L * Cdb
        Csb1L = Cgs1L * Csb
        Cgd1L = Cgs1L * Cgd
        #PHASE 2
        V2G = Vd1
        Id2 = 0.5 * uCox * W_L_2bias * (Vbias - Vss - V_th) ** 2 # Current set by Vbias
        roL2 = 1 / (lam * Id2)
        roB2 = 1 / (lam * Id2)
        ro2 = 1 / (lam * Id2)


        #Vg2 = Vov2 + Id2 * roB2 + V_th # used to find Vds2
        Vov2L = sqrt(2 * Id2 / (W_L_2L * uCox))
        #print "top transistor Vov:" + str(Vov2L)
        Vd2 = Vdd - (Vov2L + V_th)
        
        Vov2bias = sqrt(2 * Id2 / (W_L_2bias * uCox))
        Vgs2bias = Vov2bias + V_th

        Vov2 = sqrt(2 * Id2 / (W_L_2 * uCox))
        Vsb2 = 0.5 * (gam * \
            sqrt(gam ** 2 + 8 * phi + 4 * (V2G - Vov2 - V_th + gam * sqrt(2 * phi) - Vss))\
            + gam ** 2 + 2 * (V2G - Vov2 - V_th + gam * sqrt(2 * phi) - Vss))
        V_th2 = V_th + gam * (sqrt(2 * phi + Vsb2) - sqrt(2 * phi))
        Vgs2 = Vov2 + V_th2
        Vs2 = V2G - Vgs2
      
        #print "Vs2: " + str(Vs2)
        #print "Vgs2: " + str(Vgs2)
        gm2L = 2 * Id2 / (Vdd - Vd2 - V_th)
        gm2 = 2 * Id2 / Vov2
        Rl = 1/gm2L * roL2 * ro2 / (1/gm2L * roL2 + roL2 * ro2 + ro2 * 1/gm2L)
        A2 = -gm2 * Rl
        Vds2 = Vd2 - Vs2
        if Vov2 > Vds2: # in saturation
             Vbias += Vbias_inc
             continue
        # Phase 2 Capacitance Calculations:
        Cgs2 = W_L_2 * Cgate
        # Approximate other caps as ratio with Cgs1 for simplicity
        Cdb2 = Cgs2 * Cdb
        Csb2 = Cgs2 * Csb
        Cgd2 = Cgs2 * Cgd
        # Other node of stage 1/Stage2:

        Cl = Cgd1L + Cdb1L+ Cgd1 + Cdb1 + Cgs2 + Cgd2 *(1-A2)
        Ra = (RU * RD)/(RU + RD)
        Rnode2 = (roL * Ra) / (roL + Ra)
        #print Rnode2, roL, Ra
        print Cl, Rnode2
        print "HI", 1/(Cl * Rnode2 * 2 * pi * 10 ** 6) 

        # Phase 3
        V3G = Vd2

        Id3 = 0.5 * uCox * W_L_3bias * (Vbias - Vss - V_th) ** 2

        roB3 = 1 / (lam * Id3)
        ro3 = 1 / (lam * Id3)
        Vov3 = sqrt(2 * Id3 / (W_L_3 * uCox))
        Vsb3 = 0.5 * (gam * \
            sqrt(gam ** 2 + 8 * phi + 4 * (V3G - Vov3 - V_th + gam * sqrt(2 * phi) - Vss))\
            + gam ** 2 + 2 * (V3G - Vov3 - V_th + gam * sqrt(2 * phi) - Vss))
        V_th3 = V_th + gam * (sqrt(2 * phi + Vsb3) - sqrt(2 * phi))
        
        Vgs3 = Vov3 + V_th3
        Vs3 = V3G - Vgs3
        Vds3 = Vdd - Vs3
        # print "Vs3: " + str(Vs3)
        # print "Vds3: " + str(Vds3)

        gm3 = 2 * Id3 / Vov3
        gmb3 = gm3 * gam / (2 * sqrt(Vsb3 + 2 * phi ))
        Rl = roB3 * ro3 / (roB3 + ro3)
        A3 = gm3 * Rl / (1 + (gm3 + gmb3) * Rl) 
        if Vov3 > Vds3:
              Vbias += Vbias_inc
              continue
        Atot = A1 * A2 * A3
        print Atot, Amax
        if Atot > Amax:
            Amax = Atot
            #print "A", Atot, A1, A2, A3
            #print "Vbias ", Vbias, "Vbiasp ", Vbiasp
        # Phase 3 Capacitance Calculations
        Cgs3 = W_L_3 * Cgate
        Cds3 = Cgs3 * Cds
        Cgd3 = Cgs3 * Cgd
        Cds3 = Cds3 * Cds  
        Vbias += Vbias_inc
    Vbiasp += Vbias_inc
