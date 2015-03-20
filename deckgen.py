#!/usr/bin/env python3

# ========================
# QuickPIC Deck Generator
# ========================

# import sys
import numpy as _np
import jinja2 as jj


def deckgen(
        plasma, bunches, box, magic, qpic,
        restart=False, verbose_qpic=False, np=1e18, verbose=False):

    
    # =====================================
    # =====================================
    # Generate Input Deck
    # =====================================
    # =====================================
    
    # =====================================
    # Path to folder holding templates
    # =====================================
    template_path = 'templates'
    template_filename = 'rpinput'
    
    # =====================================
    # Generate jinja template object
    # =====================================
    loader        = jj.FileSystemLoader(template_path)
    env           = jj.Environment(loader=loader, trim_blocks=True)
    template      = env.get_template(template_filename)
    
    # =====================================
    # Create output file for writing
    # =====================================
    output_filename = 'rpinput.jj2'
    deckfile        = open(output_filename, 'w+')
    
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
        alpha_x0        = alpha_x0, beta_x0 = beta_x0, alpha_y0 = alpha_y0, beta_y0 = beta_y0, sig_z = sig_z,
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
        z_nstep         = z_nstep, z_prof=z_prof, z_step=z_step,
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
