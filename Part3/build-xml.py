import numpy as np

import openmc
import openmc.mgxs


###############################################################################
#                 Exporting to OpenMC materials.xml File
###############################################################################

uo2 = openmc.Material(name='uo2')
uo2.add_nuclide('U235', 0.02115, 'wo')
uo2.add_nuclide('U238', 0.86032, 'wo')
uo2.add_nuclide('O16', 0.11852, 'wo')
uo2.set_density('g/cm3', 10.3)

zirconium = openmc.Material(name='zirconium')
zirconium.add_element('Zr', 1.0)
zirconium.set_density('g/cm3', 6.55)

water = openmc.Material(name='water')
water.add_nuclide('H1', 2)
water.add_nuclide('O16', 1)
water.set_density('g/cm3', 0.701)
water.add_s_alpha_beta('c_H_in_H2O')

void = openmc.Material(name='void')
void.add_nuclide('He4', 1)
void.set_density('g/cm3', 1E-10)

pyrex = openmc.Material(name='pyrex')
pyrex.add_element('B', 4.88396e-3)
pyrex.add_element('O', 4.6624e-2)
pyrex.add_element('Al', 1.7352e-3)
pyrex.add_element('Si', 1.83512e-2)
pyrex.set_density('g/cm3', 2.26)

materials = openmc.Materials((uo2, zirconium, water, pyrex, void))
materials.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml File
###############################################################################

# FUEL PIN
pitch = 1.25984

fuel_or = openmc.ZCylinder(R=0.39218)
clad_ir = openmc.ZCylinder(R=0.40005)
clad_or = openmc.ZCylinder(R=0.45720)

fuel = openmc.Cell(1, 'fuel')
fuel.fill = uo2
fuel.region = -fuel_or

gap = openmc.Cell(2, 'air gap')
gap.fill = void
gap.region = +fuel_or & -clad_ir

clad = openmc.Cell(3, 'clad')
clad.fill = zirconium
clad.region = +clad_ir & -clad_or

moderator = openmc.Cell(4, 'moderator')
moderator.fill = water
moderator.region = +clad_or

fuel_pin = openmc.Universe()
fuel_pin.add_cells((fuel, gap, clad, moderator))

# GUIDE TUBE
clad_ir = openmc.ZCylinder(R=0.56134)
clad_or = openmc.ZCylinder(R=0.60198)

inner = openmc.Cell()
inner.fill = water
inner.region = -clad_ir

clad = openmc.Cell()
clad.fill = zirconium
clad.region = +clad_ir & -clad_or

outer = openmc.Cell()
outer.fill = water
outer.region = +clad_or

guide_tube = openmc.Universe()
guide_tube.add_cells((inner, clad, outer))

# BURNABLE POISON
radii = [0.21400, 0.23051, 0.24130, 0.42672, 0.43688, 0.48387, 0.56134, 0.60198]
cyls = [openmc.ZCylinder(R=R) for R in radii]

bp_cells = []

c = openmc.Cell()
c.region = -cyls[0]
c.fill = void
bp_cells.append(c)

mats = [zirconium, void, pyrex, void, zirconium, water, zirconium]
for i in range(len(mats)):
    c = openmc.Cell()
    c.add_surface(surface=cyls[i], halfspace=+1)
    c.add_surface(surface=cyls[i+1], halfspace=-1)
    c.fill = mats[i]
    bp_cells.append(c)

c = openmc.Cell()
c.region = +cyls[-1]
c.fill = water
bp_cells.append(c)

burn_abs = openmc.Universe()
burn_abs.add_cells(bp_cells)

# OUTSIDE
moderator = openmc.Cell()
moderator.fill = water

all_water = openmc.Universe()
all_water.add_cell(moderator)

# FUEL ASSEMBLY LATTICE
assembly = openmc.RectLattice(name='assembly')
assembly.dimension = (17, 17)
assembly.pitch = (1.26, 1.26)
assembly.lower_left = [-1.26 * 17. / 2.0] * 2

guide_tube_x = np.array([5, 8, 11, 2, 5, 8, 11, 14, 2, 5, 8,
                       11, 14, 2, 5, 8, 11, 14, 5, 8, 11])
guide_tube_y = np.array([2, 2, 2, 5, 5, 5, 5, 5, 8, 8, 8, 8,
                       8, 11, 11, 11, 11, 11, 14, 14, 14])
burn_abs_x = np.array([3, 13, 3, 13])
burn_abs_y = np.array([3, 3, 13, 13])

universes = np.empty((17, 17), dtype=openmc.Universe)
universes[:, :] = fuel_pin
universes[guide_tube_x, guide_tube_y] = guide_tube
universes[burn_abs_x, burn_abs_y] = burn_abs
assembly.universes = universes

# ROOT CELL
min_x = openmc.XPlane(x0=-10.71, boundary_type='reflective')
max_x = openmc.XPlane(x0=+10.71, boundary_type='reflective')
min_y = openmc.YPlane(y0=-10.71, boundary_type='reflective')
max_y = openmc.YPlane(y0=+10.71, boundary_type='reflective')
min_z = openmc.ZPlane(z0=-10., boundary_type='reflective')
max_z = openmc.ZPlane(z0=+10., boundary_type='reflective')

root_cell = openmc.Cell(name='root cell')
root_cell.fill = assembly
root_cell.region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z

# ROOT UNIVERSE
root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)

geometry = openmc.Geometry()
geometry.root_universe = root_universe
geometry.export_to_xml()


###############################################################################
#                   Exporting to OpenMC settings.xml File
###############################################################################

# Construct uniform initial source distribution over fissionable zones
lower_left = [-10.71, -10.71, -10]
upper_right = [10.71, 10.71, 10.]
source = openmc.source.Source(
    space=openmc.stats.Box(lower_left, upper_right, only_fissionable=True))
source.space.only_fissionable = True

# Instantiate Settings collection and export to "settings.xml"
settings = openmc.Settings()
settings.batches = 100
settings.inactive = 10
settings.particles = 10000
settings.source = source
settings.sourcepoint_write = True
settings.export_to_xml()


###############################################################################
#                     Exporting to OpenMC plots.xml File
###############################################################################

# Create a plot of the materials
plot = openmc.Plot(plot_id=1)
plot.width = [10.71 * 2] * 2
plot.basis = 'xy'
plot.color = 'mat'
plot.filename = 'materials'
plot.pixels = [2000, 2000]
plot.mask_background = [255, 255, 255]
plot.col_spec = {
    water.id:     [102, 178, 255],   # water:  blue
    zirconium.id: [ 96,  96,  96],   # zirc:   gray
    uo2.id:       [255,  75,  75],   # fuel:   red
    pyrex.id:     [  0, 153,   0],   # pyrex:  green
    void.id:      [  0,   0,   0]    # void:   white
}


# Instantiate a Materials collection and export to "materials.xml"
plots = openmc.Plots([plot])
plots.export_to_xml()


###############################################################################
#                   Create Mesh Tallies for Verification
###############################################################################

# Instantiate a tally Mesh
mesh = openmc.Mesh(name='assembly mesh')
mesh.type = 'regular'
mesh.dimension = assembly.dimension
mesh.lower_left = assembly.lower_left
mesh.width = assembly.pitch

# Instantiate mesh Filter
mesh_filter = openmc.MeshFilter(mesh)

# Fission rate mesh Tally
mesh_fiss = openmc.Tally(name='mesh fission')
mesh_fiss.filters = [mesh_filter]
mesh_fiss.scores = ['fission']
mesh_fiss.nuclides = ['U235', 'U238']

# Fine energy flux Tally
flux = openmc.Tally(name='flux')
energies = np.logspace(-2, np.log10(8e6), 1000)
flux.filters = [openmc.EnergyFilter(energies)]
flux.scores = ['flux']

# U-238 capture and fission distribcell Tally
distribcell = openmc.Tally(name='distribcell')
distribcell.filters = [openmc.DistribcellFilter(fuel.id)]
distribcell.nuclides = ['U235', 'U238']
distribcell.scores = ['absorption', 'fission']

# Resonance escape probability Tallies
therm_abs = openmc.Tally(name='thermal absorption')
therm_abs.scores = ['absorption']
therm_abs.filters = [openmc.EnergyFilter([0., 0.625])]

tot_abs = openmc.Tally(name='total absorption')
tot_abs.scores = ['absorption']

# Instantiate an 8-group structure for MGXS tallies
groups8 = openmc.mgxs.EnergyGroups()
groups8.group_edges = np.array([0.0, 0.058, 0.14, 0.28, 0.625, 4.0, 5.53e3,
                                821e3, 20e6])

# Total MGXS
tot_mgxs = openmc.mgxs.TotalXS(domain=fuel, domain_type='cell',
                               groups=groups8, by_nuclide=False)
tot_mgxs.name = 'total mgxs'

# Scattering matrix MGXS by nuclide
scatt_mgxs = openmc.mgxs.ScatterMatrixXS(
    domain=water, domain_type='cell', groups=groups8, by_nuclide=True)
scatt_mgxs.name = 'scatter matrix'

# Instantiate Tallies collection
tallies = openmc.Tallies([mesh_fiss, flux, distribcell, therm_abs, tot_abs])

# Add Tallies from each MGXS
for tally in tot_mgxs.tallies.values():
    tallies += [tally]
for tally in scatt_mgxs.tallies.values():
    tallies += [tally]

# Export Tallies collection to "tallies.xml"
tallies.export_to_xml()
