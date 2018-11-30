#!/usr/bin/python
# -*- coding: cp1252 -*-
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

# Capacitance from photodiode
Cin = 100 * 10 ** -15
# Capacitance constants (back calculation to match SPICE)
Cgate = 2.03 * 10 ** -15
# Capacitance of load
Cload = 250.e-15
# Impedance of load
Rload = 20.e3
# Ratio of parasitic capacitance
Cdb = 0.33
Csb = 0.33 # Assume Csb = Cdb 
Cgd = 0.25
FOM_Max = 0
# Sweep pvbias range:
Vbias_min = -1.9
Vbias_max = -1.4
Vbias_inc = .01
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
    vovLoad = Vdd - Vbiasp - V_th
    IdL = V_th * uCox_p * W_L_1L * vovLoad ** 2
    # RoL is given by IdL
    roL = 1/(lam * IdL)
    Vbias = Vbias_min
    print Vbiasp
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
        vov1 = -Vs1 - V_th
        gmL = (2 * IdL)/vovLoad
        gm1 = (2 * IdB1)/(vov1)
        gmb1 = gm1 * gam / (2 * sqrt(Vsb + 2 * phi))
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
        Ctot1B = Cgd1B + Cdb1B 
        Ctot1 = Cgs1 + Csb1 + Ctot1B + Cin

        BW0 = 1/(Ctot1 * R3db * 2 * pi)
        
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

        Vov2L = sqrt(2 * Id2 / (W_L_2L * uCox))
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
      
        gm2L = 2 * Id2 / (Vdd - Vd2 - V_th)
        gm2 = 2 * Id2 / Vov2
        gmb2 = gm2 * gam / (2 * sqrt(Vsb2 + 2 * phi)) # not used
        Rl = 1/gm2L * roL2 * ro2 / (1/gm2L * roL2 + roL2 * ro2 + ro2 * 1/gm2L)
        A2 = -gm2 * Rl
        Vds2 = Vd2 - Vs2
        if Vov2 > Vds2: # in saturation
             Vbias += Vbias_inc
             continue

        # Phase 2 Capacitance Calculations:
        Cgs2 = W_L_2 * Cgate
        # Approximate other caps as ratio with Cgs2 for simplicity
        Cdb2 = Cgs2 * Cdb
        Csb2 = Cgs2 * Csb
        Cgd2 = Cgs2 * Cgd

        Cgs2L = W_L_2L * Cgate
        Cdb2L = Cgs2L * Cdb
        Csb2L = Cgs2L * Csb
        Cgd2L = Cgs2L * Cgd
        # Other node of stage 1/Stage2:

        Cl = Cgd1L + Cdb1L+ Cgd1 + Cdb1 + Cgs2 + Cgd2 *(1-A2)
        Ra = (RU * RD)/(RU + RD)
        Rnode2 = Ra
        #print Rnode2, roL, Ra
      
        BW1 = 1/(Cl * Rnode2 * 2 * pi) 

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

        gm3 = 2 * Id3 / Vov3
        gmb3 = gm3 * gam / (2 * sqrt(Vsb3 + 2 * phi ))
        Rl3 = roB3 * ro3 / (roB3 + ro3)
        Rl3 = Rl3 / (Rl3 * gmb3 + 1)
        A3 = gm3 * Rl3 / (1 + gm3 * Rl3) 

        if Vov3 > Vds3:
              Vbias += Vbias_inc
              continue
        # Find final gain
        Atot = A1 * A2 * A3
        #print A1, A2, A3, Atot
        # Keep track of largest gain
        if Atot > Amax:
            Amax = Atot

        # Phase 3 Capacitance Calculations
        Cgs3 = W_L_3 * Cgate
        Cgd3 = Cgs3 * Cgd
        Csb3 = Cgs3 * Csb  

        Cgs3B = W_L_3bias * Cgate
        Cdb3B = Cgs3B * Cdb
        Cgd3B = Cgs3B * Cgd

        #Bandwidth from Stage 2 to 3
        
        C2pole = Cgd2 * (1. - 1. / A2) + Cdb2 + Cgs2L + Csb2L + Cgd3 + Cgs3 * (1 - A3) 
        R2L = roL2 / (roL2 * gm2L + 1)
        R2pole = ro2 * R2L / (ro2 + R2L)

        BW2 = 1/(C2pole * R2pole * 2 * pi)


        # Last Bandwidth Calculation


        C3pole = Cgs3*(1. - 1. / A3) + Csb3 + Cgd3B + Cdb3B + Cload
        R3pole = Rload *  Rl3 / (Rload + Rl3)

        BW3 = 1/(C3pole * R3pole * 2 * pi)
        Itot = IdL + Id2 + Id3 + (Vdd - Vd1)/RU + (Vd1 - Vss)/RD
        FOM = Atot * min(BW0, BW1, BW2, BW3)/(Itot * 5)
        if min(BW0, BW1, BW2, BW3) > FOM_Max:
            FOM_Max = min(BW0, BW1, BW2, BW3)
            print FOM, min(BW0, BW1, BW2, BW3)
            print "Params", Vbias, Vbiasp
            print "Itot", Itot
            print "Gain: ", Atot
            print "BW: ", min(BW0, BW1, BW2, BW3)/10**6
        Vbias += Vbias_inc
    Vbiasp += Vbias_inc

