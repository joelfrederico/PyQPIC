#!/usr/bin/python

#==============================================
# Deckfile Template for QuickPIC Deck Generator
#==============================================

deckstr=""
deckstr+=(\
"\
---------------QuickPIC Input------------------------- \n\
&Input_File \n\
 Version = 032011  \n\
/ \n\
 \n\
--------------Pipeline Processing--------------------- \n\
"
)
if (verbose):
    deckstr+=(\
"\
Number of stages in the pipeline \n\
------------------------------------------------------ \n\
"
)
deckstr+=(\
"\
&Pipeline \n\
 Num_Stages = "+str(N_STAGES)+"   \n\
/ \n\
 \n\
--------------Simulation System----------------------- \n\
"
)
if (verbose):
    deckstr+=(\
"\
Simulation system (in unit of micron) = BOX_X * BOX_Y \n\
 * BOX_Z \n\
Total grids = (2^INDX) * (2^INDY) * (2^INDZ) \n\
Total beam particles = NPX * NPY * NPZ \n\
------------------------------------------------------ \n\
"
)
deckstr+=(\
"\
&Simulation_Sys \n\
 Box_X= "+str(box_xy)+", Box_Y= "+str(box_xy)+", Box_Z= "+str(box_z)+", \n\
 INDX = "+str(ind_xy)+", INDY = "+str(ind_xy)+", INDZ = "+str(ind_z)+"  \n\
/ \n\
 \n\
--------------Boundary Condition---------------------- \n\
"
)
if (verbose):
    deckstr+=(\
"\
Choose between 'periodic' and 'conducting'. \n\
------------------------------------------------------ \n\
"
)
deckstr+=(\
"\
&Boundary \n\
 SBoundary = 'conducting'  \n\
/ \n\
 \n\
-------------- Beams --------------------------------- \n\
"
)
if (verbose):
    deckstr+=(\
"\
NBeams = number of beams \n\
------------------------------------------------------ \n\
"
)
deckstr+=(\
"\
&Num_Beams  \n\
 NBeams = "+str(nbeam)+"  \n\
/ \n\
 \n\
-------------Beam Parameters-------------------------- \n\
"
)
if (verbose):
    deckstr+=(\
"\
BEAM_EVOLUTION : turn on/off beam push \n\
MIN_BEAM_PARTICLE=minimal number of beam particles in  \n\
each processor. \n\
NPX, NPY, NPZ : NPX*NPY*NPZ is the total number of  \n\
particles for this beam  \n\
Charge = charge of beam particle, in unit of e. \n\
Mass = mass of beam particle, in unit of electron mass.  \n\
GAMMA = lorentz factor \n\
Num_Particle = Number of beam particle. \n\
VDX(Y&Z) = drift velocity of the beam, in unit of c \n\
Init_Routine : specify which init routine to use. \n\
             1 :  tri-gaussian random initializtion \n\
             2 :  bi-gaussian in x and y, piecewise  \n\
                  linear in z \n\
             3 :  bi-gaussian in x and y, piecewise  \n\
                  linear in z, random initialization \n\
             4 :  arbitrary 3D profile specified by  \n\
                  the BEAM_PROFILE file, parameter  \n\
                  array is ignored.  \n\
             5 : twiss parameter initialization for  \n\
                 transverse phase space, gaussian for  \n\
                 longitudinal profile. \n\
Parameter_Array = parameters for the init routine. \n\
Parameter_Array(1,:) = (Center_X,Center_Y,Center_Z)  \n\
                 = Position of the center of the beam \n\
Parameter_Array(2,:) = \n\
  Init_Routine=1 :  (Sigma_X, Sigma_Y, Sigma_Z) \n\
                        Sigma_X(Y) in micron  \n\
  Init_Routine=2,3 : (Sigma_X, Sigma_Y, Size_of_Profile_Array) \n\
                     Maximum size = 500 \n\
  Init_Routine=5 : (Alpha_X, Beta_X, Alpha_Y, Beta_Y, Sigma_Z) \n\
Parameter_Array(3,:) = \n\
  Init_Routine=1,5 : (EMITTANCE_X, EMITTANCE_Y, ENERGY_DIFF)   \n\
                   Normalized emittance of the beam in unit of  \n\
                   mm*mrad, thermal velocity of the beam =  \n\
                   emittance/(gamma*sigma) \n\
                   ENERGY_DIFF = DELTA_GAMMA/GAMMA, logitudinal  \n\
                   thermal velocity of the beam is ENERGY_DIFF/GAMMA \n\
  Init_Routine=2,3: (EMITTANCE_X, EMITTANCE_Y, ENERGY_DIFF, \n\
                   Centroid_C2X,Centroid_C1X,Centroid_C0X,Centroid_C2Y, \n\
                   Centroid_C1Y,Centroid_C0Y) \n\
                   Beam centroid is described by parabolic function \n\
                   Centroid_C2*(Z-Z0)^2+Centroid_C1*(Z-Z0)+Centroid_C0 \n\
                   Here Z and Z0 are in unit of micron, the code wil  \n\
                   convert Centroid_C2(1&0) into the unit in the \n\
                   slmulation. \n\
Parameter_Array(4,:) =  \n\
   Init_Routine=1 : (Centroid_C2X,Centroid_C1X,Centroid_C0X) \n\
   Init_Routine=2,3 : f(x(i)) in beam profile, arb. unit \n\
   Init_Routine=5 : not used  \n\
Parameter_Array(5,:) = \n\
   Init_Routine=1 : (Centroid_C2Y,Centroid_C1Y,Centroid_C0Y) \n\
   Init_Routine=2,3 : x(i) in beam profile, in micron  \n\
   Init_Routine=5 : not used  \n\
Use_Shifter = Shift particles' transverse position after  \n\
   initialization \n\
Shifter_Nsec = number of sections of the displacements \n\
Shifter_Parameter(1,:) = displacement in x direction \n\
Shifter_Parameter(2,:) = displacement in y direction \n\
Shifter_Parameter(3,:) = z position of each section \n\
Use_Destroyer = particle destroyer \n\
Destroyer_NCriteria = number of criteria \n\
Destroyer_Criteria(1,:) = dimension to operate on \n\
   (1:X, 2:Y, 3:Z, 4:Px, 5:Py, 6:Pz) \n\
Destroyer_Criteria(2,:) = lower bound \n\
Destroyer_Criteria(3,:) = upper bound \n\
   units are in micron(for XYZ) or mc(for Px,Py,Pz) \n\
   inbound particles will be destroyed! \n\
------------------------------------------------------ \n\
"
)
for i in range(0,nbeam):
    deckstr+=(\
"\
&Beam \n\
 BEAM_EVOLUTION = ."+str(BEAM_EVO)+". \n\
 MIN_BEAM_PARTICLE = 5000000  \n\
 NPX = "+str(NP_xy)+", NPY = "+str(NP_xy)+", NPZ = "+str(NP_z)+" \n\
 Charge = -1.0 \n\
 Mass = 1.0 \n\
 Gamma = "+str(gam[i])+" \n\
 Num_Particle = "+str(Q[i])+" \n\
 VDX = "+str(v_x[i])+", VDY = "+str(v_y[i])+", VDZ = "+str(v_z[i])+" \n\
 Init_Routine =5 \n\
 BEAM_PROFILE = 'test.hdf'  \n\
 QUIET_START = .true.  \n\
 Parameter_Array(1:1,1:3) = "+str(C_x[i])+","+str(C_y[i])+","+str(C_z[i])+" \n\
 Parameter_Array(2:2,1:5) = "+str(alpha_x0[i])+","+str(beta_x0[i])+","+str(alpha_y0[i])+","+str(beta_y0[i])+","+str(sig_z[i])+" \n\
 Parameter_Array(3:3,1:3) = "+str(en_x[i])+","+str(en_y[i])+","+str(dp[i])+" \n\
 Parameter_Array(4:4,1:3) = 0.,0.,0. \n\
 Parameter_Array(5:5,1:3) = 0.,0.,0. \n\
 Use_Shifter = .false.                           \n\
 Shifter_Nsec = 4   \n\
 Shifter_Parameter(1:1,1:4) = 0.,0.,0.,0.  \n\
 Shifter_Parameter(2:2,1:4) = 0.,0.,0.,0. \n\
 Shifter_Parameter(3:3,1:4) = 0.,0.,0.0,0.0  \n\
 Use_Destroyer = .true. \n\
 Destroyer_NCriteria = 5  \n\
 Destroyer_Criteria(1:1,1:5)=1,1,2,2,6 \n\
 Destroyer_Criteria(2:2,1:5)=0,"+str(box_xy-2)+",0,"+str(box_xy-2)+",0 \n\
 Destroyer_Criteria(3:3,1:5)=2,"+str(box_xy)+",2,"+str(box_xy)+",100 \n\
 Use_Radiation_Damping = .false. \n\
/ \n\
 \n\
"
)
deckstr+=(\
"\
--------------laser_input ---------------------------- \n\
&laser_input \n\
 laser_on = .false. \n\
/ \n\
 \n\
--------------plasma species------------------------ \n\
"
)
if (verbose):
    deckstr+=(\
"\
Nspecies: total number of plasma species  \n\
Plasma_Density: density for normalization, \n\
                in unit of cm-3 \n\
                not necessarily density of one species \n\
---------------------------------------------------- \n\
"
)
deckstr+=(\
"\
&Plasma \n\
 Nspecies= "+str(N_PREFORMED)+" \n\
 Nneutrals= "+str(N_NEUTRAL)+" \n\
 Plasma_Density= "+str(np)+" \n\
/ \n\
 \n\
------------Plasma Parameters------------------------- \n\
"
)
if (verbose):
    deckstr+=(\
"\
LOAD_BALANCE_TH = threshold value for load balancing. \n\
NP2: NP2*NP2 is the number of simulation particles in  \n\
    one 2D slice. \n\
Charge = charge of plasma particle, in unit of e. \n\
Mass = mass of plasma particle, in unit of electron mass. \n\
VT2X(Y) = thermal velocity of the plasma electrons, in  \n\
    unit of c \n\
Non_Neutral_Factor = - fixed ion density/electron density, \n\
    Non_Neutral_Factor = 1 for neutral plasma \n\
    Non_Neutral_Factor = 0 for pure electron cloud \n\
    Effective only when conducting boundary condition \n\
    is set. \n\
Profile_type: 0 - uniform, density = 1  \n\
                  (nomalized to the Plasma_Density) \n\
              1 - linear, density = 1+p1*(x/p2-p3) \n\
              2 - sine, density = 1+p1*sin(x*p2-p3) \n\
              3 - gaussian, density = 1+p1*exp(-((x-p2)/p3)**2)  \n\
              18 - hollow channel, density = 0 (r<p1) or p2 (r>p1) \n\
              19 - circle, density = p2 (r<p1) or 0 (r>p1)  \n\
              20 - half space, density = 0 (right) or 1 (left) \n\
              21 - piecewise, density = n(r) \n\
argx1/2/3: arguments for uniform, linear, sine, gaussian, hollow,  \n\
           circle, half space profiles. \n\
           uniform: argx1/2/3 not used. \n\
           linear:  p1=argx1, p2=argx2, p3=argx3 \n\
           sine:  p1=argx1, p2=argx2, p3=argx3 \n\
           gaussian:  p1=argx1, p2=argx2, p3=argx3 \n\
           hollow/circle: p1=argx1 (micron), p2=argx2  \n\
           half: argx1/2/3 not used            \n\
Prof_Paras_Nsec = number of points in the piecewise function. Max=100 \n\
Prof_Parameters(1,1:100): n(r) for piecewise profile, n(r>box/2) is  \n\
     forced to be 0. \n\
Prof_Parameters(2,1:100): r for piecewise profile, in micron \n\
Density_Variation: Allow density variation in s, which is  \n\
    the propagation distance. \n\
Density_Variation_Nsec: Number of sections of piece-wise  \n\
    linear function describing density variation, max=100. \n\
Density_Variation_Fs: values of piece-wise linear function.  \n\
    These are the density ratios with respect to Plasma_Density.  \n\
Density_Variation_s: corresponding propagation distances (in  \n\
    micron). \n\
------------------------------------------------------ \n\
"
)
deckstr+=(\
"\
&Species \n\
 LOAD_BALANCE_TH = -0.08 \n\
 NP2 = "+str(NP2)+"  \n\
 Charge = -1.0 \n\
 Mass = 1.0 \n\
 VT2X=0.0, VT2Y=0.0 \n\
 Non_Neutral_Factor = 1.0  \n\
 Profile_type = "+str(plas_prof)+" \n\
 argx1 = "+str(plas_p1)+" \n\
 argx2 = "+str(plas_p2)+" \n\
 argx3 = "+str(plas_p3)+" \n\
 nsrand = 0  \n\
 Prof_Nsec = 9  \n\
 Prof_Parameter(1,1:9) = 0.1,0.1,1,1,0.05,0.05,1,1,0 \n\
 Prof_Parameter(2,1:9) = 0,20,20.1,40,40.1,60,60.1,90,95 \n\
 Density_Variation=."+str(dense_var)+". \n\
 Density_Variation_Nsec="+str(z_nstep)+" \n\
 Density_Variation_Fs(1:"+str(z_nstep)+") = "+' '.join(map(str,z_prof))+" \n\
 Density_Variation_s(1:"+str(z_nstep)+") = "+' '.join(map(str,z_step))+" \n\
/ \n\
&Neutral \n\
 Neutral_gas = "+str(Z_PLAS)+" \n\
 Neutral_z = "+str(MAX_ION)+" \n\
/ \n\
 \n\
------------Simulation time--------------------------- \n\
"
)
if (verbose):
    deckstr+=(\
"\
TEND = Total time, DT = TimeStep  \n\
In unit of 1/Omega_p. \n\
------------------------------------------------------ \n\
"
)
deckstr+=(\
"\
&Simulation_time \n\
 TEND = "+str(TEND)+", DT = "+str(DT)+" \n\
/  \n\
 \n\
------------ Diagnostic ------------------------------ \n\
"
)
if (verbose):
    deckstr+=(\
"\
DFPSI, DFPHI, DFQEB, DFQEP, DFVP, DFCHI, DFJP, DFJB,  \n\
DFE, DFB  are the intervals in unit  \n\
of timestep to output PSI, PHI, beam and plasma  \n\
density, ponderomotive potential, CHI, plasma current, \n\
beam current, E field and B field respectively. \n\
DF*SLICE specify the interval for outputing 2D slices \n\
of the data. \n\
PHI(PSI,QEB,QEP)X0, if not zero, specify which Y-Z  \n\
slice to dump.  \n\
PHI(PSI,QEB,QEP)Y0, if not zero, specify which X-Z  \n\
slice to dump. \n\
PHI(PSI,QEB,QEP)Z0, if not zero, specify which X-Y  \n\
slice to dump. \n\
BC_DIAG_RES specify the number of slices along Z  \n\
direction for beam centroid calculation. \n\
------------------------------------------------------ \n\
"
)
deckstr+=(\
"\
&Potential_Diag \n\
 DFPHI=0, \n\
 DFPHISLICE=0, PHIX0=0 ,PHIY0=0, PHIZ0=0, \n\
 DFPSI=0,  \n\
 DFPSISLICE=0, PSIX0=0,PSIY0=0, PSIZ0=0 \n\
/ \n\
 \n\
&Ponderomotive_Potential_Diag \n\
 DFVP=0, \n\
 DFVPSLICE=0, VPX0=0, VPY0=0, VPZ0=0 \n\
/ \n\
 \n\
&Chi_Diag \n\
 DFCHI=0, \n\
 DFCHISLICE=0, CHIX0=0, CHIY0=0, CHIZ0=0 \n\
/ \n\
 \n\
&Current_Diag \n\
 DFJP=0,  \n\
 DFJPSLICE=0, JPX0=0, JPY0=0, JPZ0=0 \n\
 DFJB=0,  \n\
 DFJBSLICE=0, JBX0=0, JBY0=0, JBZ0=0 \n\
/ \n\
 \n\
&Field_Diag \n\
 DFE=0,  \n\
 DFESLICE= "+str(DFE)+", EX0= "+str(EX0)+", EY0= "+str(EY0)+", EZ0= "+str(EZ0)+" \n\
 DFB=0,  \n\
 DFBSLICE= "+str(DFB)+", BX0= "+str(BX0)+", BY0= "+str(BY0)+", BZ0= "+str(BZ0)+" \n\
/ \n\
 \n\
&Beam_Diag \n\
 DFQEB=0,  \n\
 DFQEBSLICE= "+str(DFQEB)+" , QEBX0= "+str(QEBX0)+", QEBY0= "+str(QEBY0)+",  QEBZ0= "+str(QEBZ0)+", \n\
 DFBC=0, BC_DIAG_RES= "+str(BC_RES)+" \n\
/ \n\
 \n\
&Plasma_Diag \n\
 DFQEP=0,  \n\
 DFQEPSLICE= "+str(DFQEP)+" , QEPX0= "+str(QEPX0)+", QEPY0= "+str(QEPY0)+",  QEPZ0= "+str(QEPZ0)+" \n\
/ \n\
 \n\
------------ Diagnostic ------------------------------ \n\
"
)
if (verbose):
    deckstr+=(\
"\
 DUMP_PHA: switch to turn on phase space diagnostics \n\
 DFPHA:  intevals in unit of timestep for dumping phase \n\
 space \n\
 DSAMPLE :  spacing of sampling \n\
------------------------------------------------------ \n\
"
)
deckstr+=(\
"\
&Beam_Phase_Space_Diag \n\
 DUMP_PHA_BEAM=."+str(DUMP_PHA_BEAM)+"., DFPHA_BEAM= "+str(DFPHA_BEAM)+",  \n\
 DSAMPLE_BEAM = "+str(DSAMPLE_BEAM)+" \n\
/ \n\
 \n\
&Plasma_Phase_Space_Diag \n\
 DUMP_PHA_PLASMA = ."+str(DUMP_PHA_PLASMA)+"., DFPHA_PLASMA = "+str(DFPHA_PLASMA)+", \n\
 DSAMPLE_PLASMA = "+str(DSAMPLE_PLASMA)+" \n\
/ \n\
 \n\
------------ Restart file ---------------------------- \n\
"
)
if (verbose):
    deckstr+=(\
"\
READ_RST_FILE specify a restart run and  RST_TIMESTEP  \n\
which timestep to begin the restart run \n\
DUMP_RST_FILE control restart file dumping and DFRST \n\
is the dumping frequency \n\
------------------------------------------------------ \n\
"
)
deckstr+=(\
"\
&Restart_File \n\
 READ_RST_FILE = ."+str(READ_RST)+"., RST_TIMESTEP = "+str(RST_START)+"  \n\
 DUMP_RST_FILE = ."+str(DUMP_RST)+".,  DFRST= "+str(DRST_STEP)+" \n\
/ \n\
 \n\
------------Optimization Coefficents---------------- \n\
"
)
if (verbose):
    deckstr+=(\
"\
INTERNAL DATA. DO NOT CHANGE! \n\
------------------------------------------------------ \n\
"
)
deckstr+=(\
"\
&Optimization \n\
 INORDER = 1, POPT = 1, DOPT = 2, DJOPT = 1  \n\
 SORTIME_2D = 25, SORTIME_3D = 25 \n\
/ \n\
 \n\
------------Debug------------------------------------- \n\
"
)
if (verbose):
    deckstr+=(\
"\
Debug options \n\
------------------------------------------------------ \n\
"
)
deckstr+=(\
"\
&debug \n\
 MAX_ITER = 2, FAC_EXY = 1., FAC_BXY = 1., FAC_AZ = 1,  \n\
 FAC_BZ = 1, C_DIF = 1 , J_DIF = 1, VERBOSE = "+str(VERBOSE)+"  \n\
/ \n\
"
)
