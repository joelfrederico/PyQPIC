#!/usr/bin/env python3

# ========================
# QuickPIC Deck Generator
# ========================
import numpy as _np
import jinja2 as jj
from . import classes as _cl


def deckgen(
        plasma, bunches, box, magic, qpic, verbose=True, filename='rpinput'
        ):
    beampos = _cl.BeamPositioning()
    phas_samp = _cl.PhaseSpaceSampling()
    # =====================================
    # =====================================
    # Generate Input Deck
    # =====================================
    # =====================================
    
    # =====================================
    # Path to folder holding templates
    # =====================================
    template_filename = 'rpinput'
    
    # =====================================
    # Generate jinja template object
    # =====================================
    loader        = jj.PackageLoader('PyQPIC', 'resources/templates')
    env           = jj.Environment(loader=loader, trim_blocks=True)
    template      = env.get_template(template_filename)
    
    # =====================================
    # Create output file for writing
    # =====================================
    deckfile        = open(filename, 'w+')
    
    # =====================================
    # Use jinja template object to create
    # a stream with substitutions in it
    # =====================================
    stream = template.stream(
        qpic              = qpic,
        magic             = magic,
        bunches           = bunches,
        plasma            = plasma,
        box               = box,
        beampos           = beampos,
        phas_samp         = phas_samp,
        verbose           = verbose,
        nbeam             = _np.size(bunches)
        )
    # =====================================
    # Write stream contents to file
    # =====================================
    stream.dump(deckfile)
    
    # =====================================
    # Close file
    # =====================================
    deckfile.close()
