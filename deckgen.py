#!/usr/bin/env python

#========================
# QuickPIC Deck Generator
#========================

import sys
from numpy import *
import jinja2 as jj

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
wp     = (5.64E4)*sqrt(np) # rad/s
# characteristic length
cwp    = c/wp                   # cm
# wavenumber
kp     = 1/cwp                  # 1/cm

#-------------------
# Calculated Peak Current
#-------------------
# peak current of drive bunch
Lambda = Q[0]*qe*c/(sqrt(2*pi)*sig_z[0]*um2cm) # A

#-------------------
# Calculated Plasma Bubble Size
#-------------------
# max bubble radius
R_bub    = 2.58*sqrt((Lambda/(qe*c))*(1/np))*cm2um # um
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
        sig_x[i] = sqrt(en_x[i]*(cm2um/kp)*sqrt(2/gam[i]))
    if (sig_y_matched[i]):
        sig_y[i] = sqrt(en_y[i]*(cm2um/kp)*sqrt(2/gam[i]))
# auto-match emittance to plasma
for i in range(0,nbeam):
    if (en_x_matched[i]):
        en_x[i] = pow(sig_x[i],2)*(kp/cm2um)*sqrt(gam[i]/2)
    if (en_y_matched[i]):
        en_y[i] = pow(sig_y[i],2)*(kp/cm2um)*sqrt(gam[i]/2)

#-------------------
# Initial beam parameters
#-------------------
for i in range(0,nbeam):
    beta_x.append   ( pow(sig_x[i],2)/en_x[i]/gam[i] )
    beta_y.append   ( pow(sig_y[i],2)/en_y[i]/gam[i] )
    beta_x0.append  ( beta_x[i] + pow(waist[i],2)/beta_x[i] )
    beta_y0.append  ( beta_y[i] + pow(waist[i],2)/beta_y[i] )
    alpha_x0.append ( waist[i]/beta_x[i] )
    alpha_y0.append ( waist[i]/beta_y[i] )
    sig_x0.append   ( sqrt(en_x[i]*beta_x0[i]) )
    sig_y0.append   ( sqrt(en_y[i]*beta_y0[i]) )

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
             7.0*max(max(sig_x0),max(sig_y0))) \
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
ind_xy   = int(max(round(log2(box_xy/d_grid)),6))
ind_z    = int(max(round(log2(box_z/d_grid)),6))
# limit number of cells in each dimension
if (ind_xy>9):
    ind_xy=9
if (ind_z>8):
    ind_z=8
# ensure ind_z > ind_xy
if (ind_z <= ind_xy):
    ind_z = int(ind_xy+1)
N_cell   = int(pow(2,ind_xy+ind_xy+ind_z))
# box volume
V_box    = box_xy*box_xy*box_z # um^3
# beam particle density (power of 2)
# default "safe" values->
NP_xy    = int(pow(2,8))
NP_z     = int(pow(2,7))
# dynamically assigned values->
#NP_xy    = int(pow(2,ind_xy+NP_xy_fact))
#NP_z     = int(pow(2,ind_z))
# plasma particle density (power of 2)
NP2      =round(sqrt(pow(pow(2,ind_xy),2)*max(np_cell,np_min)))
# beam phase space sampling particle count
ind_beam_pha = int(round(log2(samp_beam_pha_N)))
NP_beam_pha = int(pow(2,2*ind_xy+ind_z-ind_beam_pha))
# plasma phase space sampling particle count
ind_plas_pha = int(round(log2(samp_plas_pha_N)))
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
DT = round(DT_fact*sqrt(2.*gam[0])) # 1/wp
# total simulation time
if (is_multistep):
    TEND = floor(L_sim/cwp) # 1/wp
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
# initialize transverse plasma parameters
plas_p1 = 1
plas_p2 = 1
plas_p3 = 1
# set transverse profile code
# --------
#  0: flat
#  3: gaussian ramp
# 19: circle
# 21: piecewise (custom)
# --------
# flat transverse profile:
# np(r) = np0
if   (plasma_trans_geom=='flat'):
    plas_prof = int(0)
# gaussian transverse  profile:
# np(r) = np0(1+p1*exp(-((r-p2)/p3)^2))
# p1: height of peak w.r.t. np0
# p2: offset of peak
# p3: sqrt(2)*sig
elif (plasma_trans_geom=='gauss'):
    plas_prof = int(3)
    plas_p1   = -1*ramp_dir
    plas_p2   = ((1-ramp_dir)/2)*L_sim
    plas_p3   = 2*pow(ramp_width,2)
# circular filament:
# np(r) = np0*p2(r>p1) or 0(r<p1)
elif (plasma_trans_geom=='circle'):
    plas_prof = int(19)
    plas_p1 = plas_radius
    plas_p2 = 1.0

# longitudinal density profile
# --------
# initialize longitudinal density parameters
z_max   = 1
z_nstep = int(1)
z_step  = [0]
z_prof  = [0]
# --------
# define gaussian ramp function
def gaussramps(A,z1,sig1,z2,sig2,z):    
    return A*(
    (z<z1)*exp(-pow(z-z1,2)/(2*pow(sig1,2)))+
    (z>=z1)*(z<=z2)+
    (z>z2)*exp(-pow(z-z2,2)/(2*pow(sig2,2))))
# --------
# flat longitudinal profile
if   (plasma_long_geom=='flat'):
    dense_var = "false"
# gaussian ramps profile
elif (plasma_long_geom=='gauss'):
    dense_var = "true"
    z_max   = TEND*cwp*cm2um # um
    z_nstep = int(min(100,floor(TEND/DT)))
    z_step  = linspace(0,z_max,z_nstep) # um
    z_prof  = gaussramps(1,flat_start,upramp_sig,
              flat_end,dnramp_sig,z_step) # np

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

#  deckwrite=open('deckwrite.py','r')
#  exec(deckwrite)
#  deckwrite.close()

# =====================================
# Path to folder holding templates
# =====================================
template_path = 'templates'
template_filename = 'rpinput'

# =====================================
# Generate jinja template object
# =====================================
loader        = jj.FileSystemLoader(template_path)
env           = jj.Environment(loader=loader,trim_blocks=True)
template      = env.get_template(template_filename)

# =====================================
# Create output file for writing
# =====================================
output_filename = 'rpinput.jj2'
deckfile        = open(output_filename,'w+')

# =====================================
# Use jinja template object to create
# a stream with substitutions in it
# =====================================
stream = template.stream(
        verbose         = verbose,
        N_STAGES        = N_STAGES,
        box_xy          = box_xy, box_z = box_z,
        ind_xy          = ind_xy, ind_z = ind_z,
        nbeam           = nbeam,
        BEAM_EVO        = BEAM_EVO,
        NP_xy           = NP_xy, NP_z = NP_z,
        gam             = gam,
        Q               = Q,
        v_x             = v_x , v_y = v_y , v_z = v_z ,
        C_x             = C_x , C_y = C_y , C_z = C_z ,
        alpha_x0        = alpha_x0, beta_x0 = beta_x0, alpha_y0 = alpha_y0, beta_y0 = beta_y0, sig_z=sig_z,
        en_x            = en_x, en_y = en_y, dp = dp,
        N_PREFORMED     = N_PREFORMED,
        N_NEUTRAL       = N_NEUTRAL,
        np              = np,
        NP2             = NP2,
        plas_prof       = plas_prof,
        plas_p1         = plas_p1,
        plas_p2         = plas_p2,
        plas_p3         = plas_p3,
        dense_var       = dense_var,
        z_nstep         = z_nstep, z_prof=z_prof,
        Z_PLAS          = Z_PLAS,
        MAX_ION         = MAX_ION,
        TEND            = TEND, DT = DT,
        DFE             = DFE   , EX0=EX0     , EY0=EY0     , EZ0=EZ0     ,
        DFB             = DFB   , BX0=BX0     , BY0=BY0     , BZ0=BZ0     ,
        DFQEB           = DFQEB , QEBX0=QEBX0 , QEBY0=QEBY0 , QEBZ0=QEBZ0 ,
        DFQEP           = DFQEP , QEPX0=QEPX0 , QEPY0=QEPY0 , QEPZ0=QEPZ0 ,
        BC_RES          = BC_RES,
        DUMP_PHA_BEAM   = DUMP_PHA_BEAM   , DFPHA_BEAM=DFPHA_BEAM     ,
        DUMP_PHA_PLASMA = DUMP_PHA_PLASMA , DFPHA_PLASMA=DFPHA_PLASMA ,
        DSAMPLE_BEAM    = DSAMPLE_BEAM,
        DSAMPLE_PLASMA  = DSAMPLE_PLASMA,
        READ_RST        = READ_RST , RST_START=RST_START ,
        DUMP_RST        = DUMP_RST , DRST_STEP=DRST_STEP ,
        VERBOSE         = VERBOSE
        )
# =====================================
# Write stream contents to file
# =====================================
stream.dump(deckfile)

# =====================================
# Close file
# =====================================
deckfile.close()
# End Input Deck
# - - - - - - - - - - - - - - - - - - - - -

#------------------------
# Print Input Deck to File
#------------------------
# NOTE: deckfile must be called 'rpinput'
#  deckfile=open('rpinput','w')
#  deckfile.write(deckstr)
#  deckfile.close()
