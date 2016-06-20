CMD_EXAMPLE = """Example

    # open notebook and display pdb file
    nglview 1tsu.pdb
    
    # open notebook and display trajectory
    nglview 1tsu.parm7 -c traj.nc

    # open notebook and display trajectory by reading all files ending with .nc
    # make sure to use quote " "
    nglview 1tsu.parm7 -c "*.nc"
    
    # open notebook and copy myscript.py content to 1st cell
    nglview myscript.py
    
    # open my_notebook.ipynb notebook and run 1st cell
    nglview my_notebook.ipynb
    
    # running Jupyter notebook remotely
    nglview 1tsu.parm7 -c traj.nc --remote
"""
