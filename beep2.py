
#!/usr/bin/python
# -*- coding: cp1252 -*-
import random
import collections
from math import *
import sys
from math import *

verbose = False
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
Cox = 2.3 * 10 ** -15
Cov = .5 * 10 ** - 15
# Capacitance of load
Cload = 250.e-15
# Impedance of load
Rload = 20000
# Ratio of parasitic capacitance
Cdb = 0.33
Csb = 0.33 # Assume Csb = Cdb 
Cgd = 0.25
FOM_Max = 0
# Sweep pvbias range:
Vbias_min = -1.8
Vbias_max = -.1
Vbias_inc = .01
Vbias = -1.3

Vbiasp_min = .1
Vbiasp_max = 1.8
Vbiasp = 1.8#Vbiasp_min

RU = 35000.
RD = 300000.
Rinc = 1000000000.
Rmax = 200000.

W_1 = 10.
W_1L = 2.
W_1bias = 9.
L_1 = 1.
L_1L = 2.
L_1bias = 2.

W_L_1L = W_1L/L_1L
W_L_1bias = W_1bias/L_1bias
W_L_1 = W_1/L_1

W_2 = 5.
W_2L = 2.2
W_2bias = 2.
L_2 = 1.
L_2L = 2.
L_2bias = 2.4
W_L_2L = W_2L/L_2L
W_L_2bias = W_2bias/L_2bias
W_L_2 = W_2/L_2



W_3 = 28
W_3bias = 8.2
L_3 = 1.
L_3bias = 2.
W_L_3 = W_3/L_3
W_L_3bias = W_3bias/L_3bias


Amax = -1000.
# Loop through Vbiasp, determine current down
while(RU < Rmax):
    vovLoad = Vdd - Vbiasp - V_th
    IdL = .5 * uCox_p * W_L_1L * vovLoad ** 2.
    # RoL is given by IdL
    roL = 1./(lam/L_1L * IdL)
    RD = 1000.
    #Loop through bottom biasing
    while(RD < Rmax):
        vovB = Vbias - Vss - V_th    
        IdB1 = .5 * uCox * W_L_1bias * vovB ** 2.
        if verbose: print "IdB1: ", IdB1
        vov1 = sqrt(2. * IdB1 / (W_L_1 * uCox))  
        if verbose: print vov1     
        roB = 1./(lam/(L_1bias) * IdL)
        # Body Effect Time
        # This equation is a taylor expansion around vov1 = 0
        #Vsb = 1.60598 - .837037 * vov1
        #if verbose: print "Vsb: ", Vsb
        #Vs1 = Vss - Vsb
        #print Vs1
        x1 = -vov1 - V_th
        Vs1 = 2.5033e-33 * (399472499647113023775870099128320 * x1\
            - sqrt(57448180070752583817517957110915703992175254059700014444703645696 * x1\
            + 225579259037095599581293304753759236354804847248000062787317399552) + \
            286284489381129659441783532158976)
        Vsb = Vs1 - Vss
        Ires = IdL - IdB1
        # Vd1 = Current through it * R tot
        Rtot = RU * RD / (RU + RD)
        #Rtot = Rtot * roL/ (Rtot + roL)
        Vd1 = (Ires + Vdd/RU + Vss/RD)*Rtot # somewhat off 0.03V
        if verbose: print "roL:",roL
        print Vd1, Vs1
        Vds1 = Vd1 - Vs1
        if verbose: print "Vds1:", Vds1
        if vov1 > Vds1:
            print vov1, Vds1
            RD += Rinc
            continue

        gmL = (2. * IdL)/vovLoad
        gm1 = (2. * IdB1)/(vov1)
        if verbose: print "gm1: ", gm1
        gmb1 = gm1 * gam / (2. * sqrt(Vsb + 2. * phi))
        if verbose: print "gmb1: ", gmb1
        # Use vovPmos and Itop to determine W of M_L1
        Iu = (Vdd - Vd1)/RU
        Id = (Vd1 - Vss)/RD
        # Rtot = roL//Ru//Rd
        A1 = -Rtot

        # Bandwidth Calculations for Stage 1
        # Calculate Cgs1 from tech. parameters
        Cgs1 = 2. / 3. * W_1 * L_1 * Cox + Cov * W_1
        # Approximate other caps as ratio with Cgs1 for simplicity
        Cdb1 = Cgs1 * Cdb
        Csb1 = Cgs1 * Csb
        Cgd1 = Cgs1 * Cgd
        #rMiller to separate common gate ro and place in parallel
        # Miller Transform -> input = Z/(1-K), output = Z * 1/(1-1/K)
        rMiller = roB / (1. - Vd1/ Vs1)
        R3db = roB/(1. + gm1*roB)
        Cgs1B = 2./3. * W_1bias * L_1bias * Cox + Cov * W_1bias
        Cdb1B = Cgs1B * Cdb
        Cgd1B = Cgs1B * Cgd
        Ctot1B = Cgd1B + Cdb1B 
        Ctot1 = Cgs1 + Csb1 + Ctot1B + Cin
        BW0 = 1./(Ctot1 * R3db * 2. * pi)
        
        # Approximate other caps as ratio with Cgs1 for simplicity

        Cgs1L = 2./3.*W_1L * L_1L * Cox + Cov * W_1L
        Cdb1L = Cgs1L * Cdb
        Csb1L = Cgs1L * Csb
        Cgd1L = Cov * W_1L
        
        if verbose: print "Vd1: ", Vd1 # 0.03 off

        #PHASE 2
        V2G = Vd1
        Id2 = 0.5 * uCox * W_L_2bias * (Vbias - Vss - V_th) ** 2 # Current set by Vbias
        if verbose: print "Id2: ", Id2
        roL2 = 1. / (lam / L_2L * Id2)
        roB2 = 1. / (lam / L_2bias * Id2)
        ro2 = 1. / (lam / L_2 * Id2)
        Vov2L = sqrt(2. * Id2 / (W_L_2L * uCox))
        gm2L = 2. * Id2 / Vov2L
        # Taylor Expansion
        xl2 = Vdd - Vov2L - V_th
        Vs2L = 2.5033e-33 * (399472499647113023775870099128320 * xl2\
            - sqrt(57448180070752583817517957110915703992175254059700014444703645696 * xl2\
            + 225579259037095599581293304753759236354804847248000062787317399552) + \
            286284489381129659441783532158976)
        Vd2 = Vs2L
        
        if verbose: print "Vd2: ", Vd2
        Vth2L = V_th + gam * (sqrt(2. * phi + (Vd2 - Vss)) - sqrt(2. * phi))
        Vgs2L = Vov2L + Vth2L
        if verbose: print "Vgs2L:", Vgs2L
        Vds2L = Id2/gm2L
        Vov2 = sqrt(2. * Id2 / (W_L_2 * uCox))
        x2 = V2G - Vov2 - V_th
        Vs2 = 2.5033e-33 * (399472499647113023775870099128320 * x2\
            - sqrt(57448180070752583817517957110915703992175254059700014444703645696 * x2\
            + 225579259037095599581293304753759236354804847248000062787317399552) + \
            286284489381129659441783532158976)
        Vsb2 = Vs2 - Vss
        V_th2 = V_th + gam * (sqrt(2. * phi + Vsb2) - sqrt(2. * phi))
        Vgs2 = Vov2 + V_th2
        if verbose: print "Vgs2: ", Vgs2
        if verbose: print "Vs2: ", Vs2
        gm2 = 2. * Id2 / Vov2
        gmb2 = gm2 * gam / (2 * sqrt(Vsb2 + 2 * phi)) # not used
        gmb2L = gm2L * gam / (2. * sqrt(Vd2 - Vss + 2. * phi))
        Rl = 1./gm2L * roL2 * ro2 / (1./gm2L * roL2 + roL2 * ro2 + ro2 * 1./gm2L)
        Rl = Rl / (Rl * gmb2L + 1.)
        Rl = Rl / (Rl * gmb2 + 1.)
        A2 = -gm2 * Rl
        Vds2 = Vd2 - Vs2
        if verbose: print "gm2L:", gm2
        if verbose: print "Vds2:", Vds2
        if Vov2 > Vds2: # in saturation
            print "Ha"
            RD += Rinc
            continue

        # Phase 2 Capacitance Calculations:
        Cgs2 = 2./3. * W_2 * L_2 * Cox + Cov * W_2
        # Approximate other caps as ratio with Cgs2 for simplicity
        Cdb2 = Cgs2 * Cdb
        Csb2 = Cgs2 * Csb
        Cgd2 = Cgs2 * Cgd
        Cgs2L = 2./3. * W_2L * L_2L * Cox + Cov * W_2L
        Cdb2L = Cgs2L * Cdb
        Csb2L = Cgs2L * Csb
        Cgd2L = Cov * W_2L
        # Other node of stage 1/Stage2:
        Cl = Cgd1L + Cdb1L + Cgd1 + Cdb1 + Cgs2 + Cgd2 * (1. - A2)
        Ra = (RU * RD)/(RU + RD)
        Rnode2 = Ra
        #if verbose: print Rnode2, roL, Ra
      
        BW1 = 1./(Cl * Rnode2 * 2. * pi) 

        # Phase 3
        V3G = Vd2 #0.712
        if verbose: print "V3G:", V3G # 0.1V off

        Id3 = 0.5 * uCox * W_L_3bias * (Vbias - Vss - V_th) ** 2.
        if verbose: print "Id3: ", Id3

        roB3 = 1. / (lam/L_3bias * Id3)
        ro3 = 1. / (lam/L_3 * Id3)
        Vov3 = sqrt(2. * Id3 / (W_L_3 * uCox))
        if verbose: print "Vov3: ", Vov3
        x3 = V3G - Vov3 - V_th
        Vs3 = 2.5033e-33 * (399472499647113023775870099128320 * x3\
            - sqrt(57448180070752583817517957110915703992175254059700014444703645696 * x3\
            + 225579259037095599581293304753759236354804847248000062787317399552) + \
            286284489381129659441783532158976)

        Vsb3 = Vs3 - Vss
        V_th3 = V_th + gam * (sqrt(2 * phi + Vsb3) - sqrt(2 * phi))
        Vgs3 = Vov3 + V_th3
        Vds3 = Vdd - Vs3
        gm3 = 2. * Id3 / Vov3
        gmb3 = gm3 * gam / (2. * sqrt(Vsb3 + 2. * phi ))
        #Rl = roB3//Rload//1/gmb3
        Rl3 = roB3 * Rload / (roB3 + Rload)
        R13 = Rl3*ro3 /(ro3+Rl3)
        Rl3 = Rl3 / (Rl3 * gmb3 + 1.)
        #R13 = Rl3*Rload / (Rload+Rl3)
        A3 = gm3/(gm3 + 1./Rl3) 
        if verbose: print "Vsb3: ", Vsb3
        if verbose: print "Vs3:", Vs3
        if verbose: print "Vgs3:", Vgs3
        if verbose: print "Vds3:", Vds3
        if verbose: print "gm3:", gm3 
        if verbose: print "gmb3:", gmb3


        if Vov3 > Vds3:
              RD += Rinc
              continue
        # Find final gain
        Atot = A1 * A2 * A3
        #if verbose: print A1, A2, A3, Atot
        # Keep track of largest gain
        if Atot > Amax:
            Amax = Atot

        # Phase 3 Capacitance Calculations
        Cgs3 = 2. / 3. * W_3 * L_3 * Cox + Cov * W_3
        Cgd3 = Cov * W_3
        Csb3 = Cgs3 * Csb  
        Cgs3B = 2. / 3. * W_3bias * L_3bias * Cox + Cov * W_3bias
        Cdb3B = Cgs3B * Cdb
        Cgd3B = Cov * W_3bias
        #Bandwidth from Stage 2 to 3
       
        C2pole = Cgd2 * (1. - 1 / A2) + Cdb2 + Cgs2L + Csb2L + Cgd3 + Cgs3 * (1. - A3) 
        R2L = roL2 / (roL2 * gm2L + 1.)
        R2pole = ro2 * R2L / (ro2 + R2L)
        BW2 = 1./(C2pole * R2pole * 2. * pi)

        
        # Last Bandwidth Calculation
        C3pole = Csb3 + Cgd3B + Cdb3B + Cload + Cgs3 * (1. - 1. / A3) 
        BW3 = 1./(C3pole * Rl3 * 2. * pi)
        BW = min(BW0, BW1, BW2, BW3)
        Itot = IdL + Id2 + Id3 + (Vdd - Vd1)/RU + (Vd1 - Vss)/RD
        FOM = Atot * BW/Itot
        if FOM > FOM_Max:
            #if abs(Vs3) < .52:
                if BW > 70e6:
            #        if Atot > 35000:
                        FOM_Max = FOM
                        Params = [RU, RD]
                        BW = BW / 1e6
                        #print "Params", Vbias, Vbiasp
                        #print "Itot", Itot
                        #print "Gain: ", A1, A2, A3, Atot
                        #print "BW: ", BW/10**6
                        #print "Vout: ", Vs3
        RD += Rinc
    RU += Rinc    
print Params, BW
