#!/usr/bin/python

#========================
# QuickPIC Deck Generator
#========================

import sys
import math
import array

# - - - - - - - - - - - - - - - - - - - - -
# Read in Input Values

# take deck input file from command line
if (len(sys.argv)>1):
    inputname = sys.argv[1]
else:
    inputname = 'deckinput.py'

# load deck input file
print 'loading deck input file '+str(inputname)
deckinput=open(inputname,'r')
exec(deckinput)
deckinput.close()

# End Inputs
# - - - - - - - - - - - - - - - - - - - - -

#-------------------
# File Restart
#-------------------
# read the restart file
# if this is a restart run
if (is_restart):
    READ_RST="true"
else:
    READ_RST="false"
# dump restart files
DUMP_RST="true"

# - - - - - - - - - - - - - - - - - - - - -
# Set Simulation Parameters

#-------------------
# QPIC Debugger Verbosity
#-------------------
VERBOSE = qpic_debug

#-------------------
# Fundamental Constants
#-------------------
# speed of light
c      = 3.E10    # cm/s
# charge of electron
qe     = 1.6E-19  # C
# pi = 3.14...
pi     = math.pi
# microns to centimeters
um2cm  = 1.E-4    # cm/um
# centimeters to microns
cm2um  = 1.E4     # um/cm
# centimeters to meters
cm2m   = 0.01     # m/cm
# meters to centimeters
m2cm   = 100.     # cm/m

#-------------------
# Calculated Plasma Parameters
#-------------------
# plasma frequency
wp     = (5.64E4)*math.sqrt(np) # rad/s
# characteristic length
cwp    = c/wp                   # cm
# wavenumber
kp     = 1/cwp                  # 1/cm

#-------------------
# Calculated Peak Current
#-------------------
# peak current of drive bunch
Lambda = Q[0]*qe*c/(math.sqrt(2*pi)*sig_z[0]*um2cm) # A

#-------------------
# Calculated Plasma Bubble Size
#-------------------
# max bubble radius
R_bub    = 2.58*math.sqrt((Lambda/(qe*c))*(1/np))*cm2um # um
# bubble length
L_bub    = 2*pi*cwp*cm2um      # um
# outer radius of bubble sheath
r_sheath = 1.287*R_bub         # um

#-------------------
# Auto-Matching to Plasma
#-------------------
# auto-match spot size to plasma
for i in range(0,nbeam):
    if (sig_x_matched[i]):
        sig_x[i] = math.sqrt(en_x[i]*(cm2um/kp)*math.sqrt(2/gam[i]))
    if (sig_y_matched[i]):
        sig_y[i] = math.sqrt(en_y[i]*(cm2um/kp)*math.sqrt(2/gam[i]))
# auto-match emittance to plasma
for i in range(0,nbeam):
    if (en_x_matched[i]):
        en_x[i] = pow(sig_x[i],2)*(kp/cm2um)*math.sqrt(gam[i]/2)
    if (en_y_matched[i]):
        en_y[i] = pow(sig_y[i],2)*(kp/cm2um)*math.sqrt(gam[i]/2)

#-------------------
# "Magic" Simulation Factors
#-------------------
box_xy_fact = 1.0         # should be ~1-2
box_z_fact  = 1.0         # should be ~1-2
d_grid_fact = 0.05        # exact
NP_xy_fact  = -1          # exact
np_min      = 4           # exact
DT_fact     = (2./3.)/10. # exact

#-------------------
# Box Size and Particle Densities
#-------------------
# box size in x and y:
# use specified size if given
if (box_width>0):
    box_xy = int(round(box_width/5.)*5)      # um
# otherwise use default value
else:
    box_xy = int(round(                      \
             box_xy_fact*max(4.0*R_bub,      \
             5.0*max(max(sig_x),max(sig_y))) \
             /5.)*5)                         # um
# ensure box_xy is odd:
if not(box_xy%2):
    box_xy += 1
# box size in z
# use specified size if given
if (box_length>0):
    box_z = int(round(box_length/5.)*5)      #um
# otherwise use default value
else:
    box_z = int(round(                       \
            box_z_fact*max(2.5*L_bub,        \
            6.0*max(sig_z))                  \
            /5.)*5)                          # um
# ensure box_z is odd:
if not(box_z%2):
    box_z  += 1
# grid spacing
d_grid   = d_grid_fact*cwp*cm2um             # um
# number of cells (power of 2)
ind_xy   = int(max(round(math.log(box_xy/d_grid,2)),6))
ind_z    = int(max(round(math.log(box_z/d_grid,2)),6))
# limit number of cells in each dimension
if (ind_xy>8):
    ind_xy=8
if (ind_z>9):
    ind_z=9
# ensure ind_z > ind_xy
if (ind_z <= ind_xy):
    ind_z = int(ind_xy+1)
N_cell   = int(pow(2,ind_xy+ind_xy+ind_z))
# box volume
V_box    = box_xy*box_xy*box_z # um^3
# beam particle density (power of 2)
# default "safe" values->
NP_xy    = int(pow(2,7))
NP_z     = int(pow(2,8))
# dynamically assigned values->
#NP_xy    = int(pow(2,ind_xy+NP_xy_fact))
#NP_z     = int(pow(2,ind_z))
# plasma particle density (power of 2)
NP2      =round(math.sqrt(pow(pow(2,ind_xy),2)*max(np_cell,np_min)))
# beam phase space sampling particle count
ind_beam_pha = int(round(math.log(samp_beam_pha_N,2)))
NP_beam_pha = int(pow(2,2*ind_xy+ind_z-ind_beam_pha))
# plasma phase space sampling particle count
ind_plas_pha = int(round(math.log(samp_plas_pha_N,2)))
NP_plas_pha = int(pow(2,2*ind_xy+ind_z-ind_plas_pha))

#-------------------
# Beam Positioning
#-------------------
C_x=[]
C_y=[]
C_z=[]
# beam centered in x & y
# beam toward front of box in z
for i in range(0,nbeam):
    C_x.append(0.5*box_xy   + off_x[i]); # um
    C_y.append(0.5*box_xy   + off_y[i]); # um
    C_z.append(4.0*sig_z[0] + off_z[i]); # um

#-------------------
# Simulation Time Scales
#-------------------
# time step size
DT = round(DT_fact*math.sqrt(2.*gam[0])) # 1/wp
# total simulation time
if (is_multistep):
    TEND = math.floor(L_sim/cwp) # 1/wp
else:
    TEND = DT*1.01               # 1/wp

#-------------------
# Beam Evolution
#-------------------
if (is_multistep):
    BEAM_EVO = "true"
else:
    BEAM_EVO = "false"

#-------------------
# Plasma Species
#-------------------
if (preformed):
    N_PREFORMED = int(1)
    N_NEUTRAL   = int(0)
else:
    N_PREFORMED = int(0)
    N_NEUTRAL   = int(1)
# plasma species only used
# if plasma is not preformed
Z_PLAS  = int(plasma_z)
MAX_ION = int(max_ion_lv)

#-------------------
# Plasma Geometry
#-------------------
# initialize plasma parameters
plas_p1 = 1
plas_p2 = 1
plas_p3 = 1
# set profile code
# --------
#  0: flat
#  3: gaussian ramp
# 19: circle
# --------
# flat density profile:
if   (plasma_geom=='flat'):
    plas_prof = int(0)
# gaussian ramp:
# np(z) = np0(1+p1*exp(-((z-p2)/p3)^2))
# p1: height of peak w.r.t. np0
# p2: offset of peak
# p3: 2*sig^2
elif (plasma_geom=='gauss'):
    plas_prof = int(3)
    plas_p1   = -1*ramp_dir
    plas_p2   = ((1-ramp_dir)/2)*L_sim
    plas_p3   = 2*pow(ramp_width,2)
# circular filament:
# np(z) = np0*p2(r>p1) or 0(r<p1)
elif (plasma_geom=='circle'):
    plas_prof = int(19)
    plas_p1 = plas_radius
    plas_p2 = 1.0

#-------------------
# Simulation Output Parameters
#-------------------
# 3-D Sampling
#-------------------
# 3-D E-field sampling?
if (samp_E_3D):
    E3D  = int(1) # 1:true
else:
    E3D  = int(0) # 0:false
# 3-D B-field sampling?
if (samp_B_3D):
    B3D  = int(1) # 1:true
else:
    B3D  = int(0) # 0:false
# 3-D Beam sampling?
if (samp_beam_3D):
    QEB3D = int(1) # 1:true
else:
    QEB3D = int(0) # 0:false
# 3-D Plasma sampling?
if (samp_plas_3D):
    QEP3D = int(1) # 1:true
else:
    QEP3D = int(0) # 0:false
#-------------------
# Phase Space Sampling
#-------------------
# Beam Phase Space Sampling?
if (samp_beam_pha):
    DUMP_PHA_BEAM = "true"
else:
    DUMP_PHA_BEAM = "false"
# Plasma Phase Space Sampling?
if (samp_plas_pha):
    DUMP_PHA_PLASMA = "true"
else:
    DUMP_PHA_PLASMA = "false"
#-------------------
# Sampling Plane(s)
#-------------------
# E-field sampling plane(s)
if (samp_E_x0):
    EX0  = box_xy/2.  # um
else:
    EX0  = 0.         # none
if (samp_E_y0):
    EY0  = box_xy/2.  # um
else:
    EY0  = 0.         # none
if (samp_E_z0):
    EZ0  = box_z/2.   # um
else:
    EZ0  = 0.         # none
# B-field sampling plane(s)
if (samp_B_x0):
    BX0  = box_xy/2.  # um
else:
    BX0  = 0.         # none
if (samp_B_y0):
    BY0  = box_xy/2.  # um
else:
    BY0  = 0.         # none
if (samp_B_z0):
    BZ0  = box_z/2.   # um
else:
    BZ0  = 0.         # none
# Beam sampling plane(s)
if (samp_beam_x0):
    QEBX0 = box_xy/2. # um
else:
    QEBX0 = 0.        # none
if (samp_beam_y0):
    QEBY0 = box_xy/2. # um
else:
    QEBY0 = 0.        # none
if (samp_beam_z0):
    QEBZ0 = box_z/2.  # um
else:
    QEBZ0 = 0.        # none
# Plasma sampling plane(s)
if (samp_plas_x0):
    QEPX0 = box_xy/2. # um
else:
    QEPX0 = 0.        # none
if (samp_plas_y0):
    QEPY0 = box_xy/2. # um
else:
    QEPY0 = 0.        # none
if (samp_plas_z0):
    QEPZ0 = box_z/2.  # um
else:
    QEPZ0 = 0.        # none
#-------------------
# Sampling Periods in Time
#-------------------
if (is_multistep):
    # E-field
    DFE      = int(samp_E_dt)    # DT
    # B-field
    DFB      = int(samp_B_dt)    # DT
    # Beam
    DFQEB    = int(samp_beam_dt) # DT
    # Plasma
    DFQEP    = int(samp_plas_dt) # DT
else:
    # E-field
    DFE      = int(1)        # DT
    # B-field
    DFB      = int(1)        # DT
    # Beam
    DFQEB    = int(1)        # DT
    # Plasma
    DFQEP    = int(1)        # DT
# Beam Phase Space
if (samp_beam_pha):
    DFPHA_BEAM   = int(samp_beam_pha_dt) # DT
else:
    DFPHA_BEAM   = int(1) # DT
# Plasma Phase Space
if (samp_plas_pha):
    DFPHA_PLASMA = int(samp_plas_pha_dt) # DT
else:
    DFPHA_PLASMA = int(1) # DT

#-------------------
# Phase Space Sampling Sizes & Resolution
#-------------------
# Beam Phase Space
DSAMPLE_BEAM   = NP_beam_pha # particles
# Plasma Phase Space
DSAMPLE_PLASMA = NP_plas_pha # particles
# Beam Centroid Resolution
BC_RES         = int(beam_cent_res_dz) # z-slices

#-------------------
# Number of Simulation Stages
#-------------------
# check that total proc. < number of cells in z
#if (TOT_PROC > pow(2,ind_z)):
#    TOT_PROC = pow(2,ind_z)
# number of proc. per stage
STAGE_PROC   = max(pow(2,ind_xy-3),1)
# number of stages = total proc./proc. per stage
N_STAGES     = max(int(TOT_PROC/STAGE_PROC),1)

# End Simulation Parameters
# - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - -
# Generate Input Deck

deckwrite=open('deckwrite.py','r')
exec(deckwrite)
deckwrite.close()

# End Input Deck
# - - - - - - - - - - - - - - - - - - - - -

#------------------------
# Print Input Deck to File
#------------------------
# NOTE: deckfile must be called 'rpinput'
deckfile=open('rpinput','w')
deckfile.write(deckstr)
deckfile.close()
