from matplotlib import rcParams

# ---------------------------------------------------------------------------------------------
# Project figure style. Font sizes were doubled on 2026-09-01: the notebooks had been overriding
# the declared sizes downward to fit dense multi-panel figures, so text rendered at ~10 pt and was
# unreadable in saved PNGs. These are now the single source of truth -- do not override them
# downward in a notebook; enlarge `figsize` instead so the panels have room.
# ---------------------------------------------------------------------------------------------

rcParams['font.family']         = 'DejaVu Sans'
rcParams['font.weight']         = 'bold'  
rcParams['font.size']           = 20
# rcParams['figure.figsize']      = [15, 8]
rcParams['figure.titlesize']    = 26
rcParams['figure.titleweight']  = 'bold'
rcParams['axes.titlesize']      = 22
rcParams['axes.labelsize']      = 20
rcParams['axes.titleweight']    = 'bold' 
rcParams['axes.labelweight']    = 'bold'
rcParams['axes.linewidth']      = 2
rcParams['xtick.top']           = True
rcParams['ytick.right']         = True
rcParams['xtick.minor.visible'] = True
rcParams['ytick.minor.visible'] = True
rcParams['xtick.labelsize']     = 18  
rcParams['ytick.labelsize']     = 18  
rcParams['xtick.direction']     = 'inout'
rcParams['ytick.direction']     = 'inout'
rcParams['xtick.major.width']   = 2
rcParams['ytick.major.width']   = 2
rcParams['xtick.minor.width']   = 1.5
rcParams['ytick.minor.width']   = 1.5
rcParams['xtick.major.size']    = 6
rcParams['ytick.major.size']    = 6
rcParams['xtick.minor.size']    = 4
rcParams['ytick.minor.size']    = 4
rcParams['legend.fontsize']     = 16
rcParams['savefig.bbox']        = 'tight'
rcParams['savefig.dpi']         = 200