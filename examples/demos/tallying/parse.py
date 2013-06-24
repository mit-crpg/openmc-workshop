#!/usr/bin/env python

import os

import silomesh as sm
import statepoint

def main():

    print "Loading statepoint..."
    sp = statepoint.StatePoint('statepoint.300.binary')

    print "Parsing statepoint..."
    sp.read_results()

    print "\nAvailable Tallies:"
    for tally in sp.tallies:
        print "\tTally {}".format(tally.id)
        print "\t\tscores: {}".format(" ".join(tally.scores))
        print "\t\tfilters: {}".format(" ".join(tally.filters.keys()))
    

    ## Loop through the mesh for tally 1 and create plots
    print "Extracting Mesh Values from Tally 1..."

    tallyid = 0     # tally 1
    score = 0       # flux is zero see (order of tally.filters.keys())

    # get mesh dimensions
    meshid = sp.tallies[tallyid].filters['mesh'].bins[0]
    for i,m in enumerate(sp.meshes):
        if m.id == meshid:
          mesh = m
    nx,ny,nz = mesh.dimension

    # loop through mesh and extract values
    thermal = {}
    fast = {}
    for x in range(1,nx+1):
        for y in range(1,ny+1):
            for z in range(1,nz+1):
                val,err = sp.get_value(tallyid,
                                       [('mesh',(x,y,z)),('energyin',0)],
                                       score)
                thermal[(x,y,z)] = val
                val,err = sp.get_value(tallyid,
                                       [('mesh',(x,y,z)),('energyin',1)],
                                       score)
                fast[(x,y,z)] = val

    # write gnuplot datafile with axially-integrated values
    print "Making gnuplots ..."
    with open('meshdata.dat','w') as fh:
        for x in range(1,nx+1):
            for y in range(1,ny+1):
                thermalval = 0.
                fastval = 0.
                for z in range(1,nz+1):
                  thermalval += thermal[(x,y,z)]
                  fastval += fast[(x,y,z)]
                fh.write("{} {} {} {}\n".format(x,y,thermalval,fastval))
                
    # write gnuplot file and make plots
    with open('tmp.gnuplot','w') as fh:
      fh.write(r"""set terminal png size 800,400
set output 'fluxplot.png'
set nokey
set autoscale fix
set multiplot layout 1,2 title "Pin Mesh Axially-Integrated Flux Tally"
set title "Thermal"
plot 'meshdata.dat' using 1:2:3 with image
set title "Fast"
plot 'meshdata.dat' using 1:2:4 with image
""")
    os.system("gnuplot < tmp.gnuplot")
  
    # make 3d silo file
    print "Making 3D SILO File ..."
    sm.init_silo("fluxtally.silo")
    sm.init_mesh('tally_mesh',nx,ny,nz,
                          -85.8774, -85.8774, -81.662, # lower left of mesh
                          1.63576,1.63576,8.1662)      # width of mesh cells
    sm.init_var('flux_tally_thermal')
    for x in range(1,nx+1):
      for y in range(1,ny+1):
          for z in range(1,nz+1):
            sm.set_value(float(thermal[(x,y,z)]),x,y,z)
    sm.finalize_var()
    sm.init_var('flux_tally_fast')
    for x in range(1,nx+1):
      for y in range(1,ny+1):
          for z in range(1,nz+1):
              sm.set_value(float(fast[(x,y,z)]),x,y,z)
    sm.finalize_var()
    sm.finalize_mesh()
    sm.finalize_silo()
  
  
if __name__ == "__main__":
    main()


