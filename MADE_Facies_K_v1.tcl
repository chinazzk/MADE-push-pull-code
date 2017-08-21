# This run is for the Macrodispersion Experiment (MADE) site push pull test, and represents
# the fBm K-field under variable saturation conditions.

# The original test was conducted by Lui et al. 2010, this study uses the same material 
# properties and modeling configuration.

# All of the units for this solution are in meters.

# This is a 3-D test run

##################################################
#Import the ParFlow TCL package
##################################################
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*

pfset FileVersion 4

pfset Process.Topology.P 2
pfset Process.Topology.Q 2
pfset Process.Topology.R 1

##################################################
# Computational Grid
##################################################

# Defining the computational grid as a 50 X 50 meter grid to a
# depth of 10.75 meters. The computational grid for this simulation 
# must be the same as the domain geometery.   

pfset ComputationalGrid.Lower.X         0
pfset ComputationalGrid.Lower.Y         0
pfset ComputationalGrid.Lower.Z         0

pfset ComputationalGrid.NX              128
pfset ComputationalGrid.NY              128
pfset ComputationalGrid.NZ              640

pfset ComputationalGrid.DX              0.390625
pfset ComputationalGrid.DY              0.390625
pfset ComputationalGrid.DZ              0.016796875

##################################################
# Names of the GeomInputs
##################################################
# Defining how the model will be filled.  For the homogenous run
# there is only one value.  

pfset GeomInput.Names                   "domaininput"
pfset GeomInput.domaininput.GeomName    domain
pfset GeomInput.domaininput.InputType   Box

##################################################
# Domain Geometry
##################################################
# Because we are using a box for the domain the domain geometery
# is directly related to the domain grid. 

pfset Geom.domain.Lower.X     0
pfset Geom.domain.Lower.Y     0
pfset Geom.domain.Lower.Z     0

pfset Geom.domain.Upper.X    50
pfset Geom.domain.Upper.Y    50
pfset Geom.domain.Upper.Z    10.75

pfset Geom.domain.Patches "left right front back bottom top"

##################################################
# Permeability
##################################################
# The hydraulic conductivity values are in the indicator file "Facies_K.pfb", which is 
# a result of the fBm K field applied to the  direct push data collected by 
# Dogan et al. 2011.

pfset Geom.Perm.Names                 "domain"
pfset Geom.domain.Perm.FileName       "madepermFbm.pfb"
pfset Geom.domain.Perm.Type	      PFBFile

pfset Perm.TensorType                 TensorByGeom
pfset Geom.Perm.TensorByGeom.Names    "domain"
pfset Geom.domain.Perm.TensorValX     1.0d0
pfset Geom.domain.Perm.TensorValY     1.0d0
pfset Geom.domain.Perm.TensorValZ     1.0d0

##################################################
# Specific Storage
##################################################

# The value used by Lui et al., 2010 for the specific storage is 1e-5 1/m

pfset SpecificStorage.Type            Constant
pfset SpecificStorage.GeomNames       "domain"
pfset Geom.domain.SpecificStorage.Value 1.00e-6

##################################################
# Phases
##################################################
# These values are set to 1 in order to use hydraulic conductivity
# instead of permeability.

pfset Phase.Names                   "water"
pfset Phase.water.Density.Type      Constant
pfset Phase.water.Density.Value     1.0

pfset Phase.water.Viscosity.Type    Constant
pfset Phase.water.Viscosity.Value   1.0

pfset Phase.water.Mobility.Type     Constant
pfset Phase.water.Mobility.Value    1.0

pfset PhaseSources.water.Type                  Constant
pfset PhaseSources.water.GeomNames             domain
pfset PhaseSources.water.Geom.domain.Value     0.0

##################################################
# Contaminant
##################################################
pfset Contaminants.Names             ""

##################################################
# Retardation
##################################################
pfset Geom.Retardation.GeomNames      ""

###################################################
# Gravity
##################################################
pfset Gravity                         1.0

#################################################
# Timing Information
##################################################
# The total time for all of the phases adds to 460 hours.

pfset TimingInfo.BaseUnit               1.0
pfset TimingInfo.StartCount             0.0
pfset TimingInfo.StartTime              0.0
pfset TimingInfo.StopTime               460
pfset TimingInfo.DumpInterval           1
pfset TimeStep.Type			Constant
pfset TimeStep.Value			1.0

###################################################
# Porosity
##################################################
# This value has been commonly used at the MADE site for solute transport modeling.  The original value
# is from Boggs et al. 1992.  

pfset Geom.Porosity.GeomNames        "domain"

pfset Geom.domain.Porosity.Type      Constant
pfset Geom.domain.Porosity.Value     0.32

###################################################
# Domain
##################################################
pfset Domain.GeomName                domain

###################################################
# Relative Permeability
##################################################
# hydraulic properties were obtained from Schapp et al., 2001.  
# Importantly this is for soil material rather than aquifer material; however, 
# the values for material are low enough that the material is not as diffuse and this is an exceptable solution.

pfset Phase.RelPerm.Type             VanGenuchten
pfset Phase.RelPerm.GeomNames        "domain"     
  

pfset Geom.domain.RelPerm.Alpha      12.4
pfset Geom.domain.RelPerm.N          2.28
# alpha and n values are representative of sand

###################################################
# Saturation
##################################################
pfset Phase.Saturation.Type                     VanGenuchten
pfset Phase.Saturation.GeomNames                "domain"

pfset Geom.domain.Saturation.Alpha              12.4         
pfset Geom.domain.Saturation.N                  2.28
pfset Geom.domain.Saturation.SRes               0.057
pfset Geom.domain.Saturation.SSat               1.0
# alpha, n, Sres, and Ssat values are representative of sand 

###################################################
# Wells
##################################################
# The injection rate is 8.18 m^3/day (0.341 m^3/hour), 
# and the extraction rate is 7.90 m^3/day (0.329 m^3/hour)
# Importantly the well is situated in the middle of the domain.
# The well was screened from ~2.4-8.0/8.5 m below the land's surface
# This meant that from the base of the domain the well starts at 2.4 m
# and the top of the well would be at 8.5 m.  However, because the top 
# of the well screen is above the water table we decreased the total  
# length of the screen by 0.5 m for the wait and extraction phases.                 
pfset Wells.Names          "well1 well2"

pfset Wells.well1.InputType     Vertical
pfset Wells.well1.Action        Injection
pfset Wells.well1.Type          Flux
pfset Wells.well1.X             25.1953125
pfset Wells.well1.Y             25.1953125
pfset Wells.well1.ZUpper        8.0
pfset Wells.well1.ZLower        2.4
pfset Wells.well1.Method        Weighted
pfset Wells.well1.Cycle         "onoff"
pfset Wells.well1.injection.Flux.water.Value 0.313
pfset Wells.well1.nothing.Flux.water.Value 0
pfset Wells.well1.extraction.Flux.water.Value -0.329
pfset Wells.well1.injection.Saturation.water.Value 1.0
pfset Wells.well1.nothing.Saturation.water.Value 1.0
pfset Wells.well1.extraction.Saturation.water.Value 1.0 

pfset Wells.well2.InputType     Vertical
pfset Wells.well2.Action        Injection
pfset Wells.well2.Type          Flux
pfset Wells.well2.X             25.1953125
pfset Wells.well2.Y             25.1953125
pfset Wells.well2.ZUpper        8.5
pfset Wells.well2.ZLower        8.0
pfset Wells.well2.Method        Weighted
pfset Wells.well2.Cycle         "onoff"
pfset Wells.well2.injection.Flux.water.Value 0.028
pfset Wells.well2.nothing.Flux.water.Value 0
pfset Wells.well2.extraction.Flux.water.Value 0
pfset Wells.well2.injection.Saturation.water.Value 1.0
pfset Wells.well2.nothing.Saturation.water.Value 1.0
pfset Wells.well2.extraction.Saturation.water.Value 1.0 


##############################################
# Timing Cycles
##################################################
# This is another standard time cycle for steady state flow.

pfset Cycle.Names          "onoff"

pfset Cycle.onoff.Names "injection nothing extraction"
pfset Cycle.onoff.injection.Length  31
pfset Cycle.onoff.nothing.Length 19
pfset Cycle.onoff.extraction.Length   410
pfset Cycle.onoff.Repeat 1

##################################################
# Boundary Conditions
##################################################
# For the purpose of this model we have defined the edges of the 
# domain as constant head boundaries.  While the top and bottom
# of the domain are no flow (i.e. constant flux of 0) 

pfset BCPressure.PatchNames         [pfget Geom.domain.Patches]
# Geom.domain.Patches "left right front back bottom top"

pfset Patch.top.BCPressure.Type            FluxConst
pfset Patch.top.BCPressure.Cycle           "onoff"
pfset Patch.top.BCPressure.injection.Value     0.0
pfset Patch.top.BCPressure.nothing.Value       0.0
pfset Patch.top.BCPressure.extraction.Value    0.0

pfset Patch.bottom.BCPressure.Type         FluxConst
pfset Patch.bottom.BCPressure.Cycle        "onoff"
pfset Patch.bottom.BCPressure.injection.Value  0.0
pfset Patch.bottom.BCPressure.nothing.Value    0.0
pfset Patch.bottom.BCPressure.extraction.Value 0.0

pfset Patch.left.BCPressure.Type           DirEquilRefPatch
pfset Patch.left.BCPressure.Cycle          "onoff"
pfset Patch.left.BCPressure.RefGeom         domain
pfset Patch.left.BCPressure.RefPatch        bottom
pfset Patch.left.BCPressure.injection.Value   8.6
pfset Patch.left.BCPressure.nothing.Value     8.6
pfset Patch.left.BCPressure.extraction.Value  8.6

pfset Patch.right.BCPressure.Type           DirEquilRefPatch
pfset Patch.right.BCPressure.Cycle          "onoff"
pfset Patch.right.BCPressure.RefGeom         domain
pfset Patch.right.BCPressure.RefPatch        bottom
pfset Patch.right.BCPressure.injection.Value   8.6
pfset Patch.right.BCPressure.nothing.Value     8.6
pfset Patch.right.BCPressure.extraction.Value  8.6

pfset Patch.front.BCPressure.Type           DirEquilRefPatch
pfset Patch.front.BCPressure.Cycle          "onoff"
pfset Patch.front.BCPressure.RefGeom         domain
pfset Patch.front.BCPressure.RefPatch        bottom
pfset Patch.front.BCPressure.injection.Value   8.6
pfset Patch.front.BCPressure.nothing.Value     8.6
pfset Patch.front.BCPressure.extraction.Value  8.6

pfset Patch.back.BCPressure.Type           DirEquilRefPatch
pfset Patch.back.BCPressure.Cycle          "onoff"
pfset Patch.back.BCPressure.RefGeom         domain
pfset Patch.back.BCPressure.RefPatch        bottom
pfset Patch.back.BCPressure.injection.Value   8.6
pfset Patch.back.BCPressure.nothing.Value     8.6
pfset Patch.back.BCPressure.extraction.Value  8.6

##################################################
# Topography
##################################################
pfset TopoSlopesX.Type              "Constant"
pfset TopoSlopesX.GeomNames         "domain"
pfset TopoSlopesX.Geom.domain.Value 0.0

pfset TopoSlopesY.Type              "Constant"
pfset TopoSlopesY.GeomNames         "domain"
pfset TopoSlopesY.Geom.domain.Value 0.0

##################################################
# Mannings Coeffiecient
##################################################
pfset Mannings.Type                 "Constant"
pfset Mannings.GeomNames            "domain"
pfset Mannings.Geom.domain.Value    0.0

##################################################
# Initial Conditions
##################################################
# The initial conditions use the height of the water table 

pfset ICPressure.Type                 HydroStaticPatch
pfset ICPressure.GeomNames            domain
pfset Geom.domain.ICPressure.Value    8.6
pfset Geom.domain.ICPressure.RefGeom  domain
pfset Geom.domain.ICPressure.RefPatch bottom

#######################################################
# Exact solution specification for error calculations
#######################################################
pfset KnownSolution                 NoKnownSolution

##################################################
# Solver
##################################################
pfset Solver                                         Richards
pfset Solver.MaxIter                                 50000

pfset Solver.Nonlinear.MaxIter                           50
pfset Solver.Nonlinear.ResidualTol                       1e-4
pfset Solver.Nonlinear.EtaChoice                         EtaConstant
pfset Solver.Nonlinear.EtaValue                          1e-5
pfset Solver.Nonlinear.UseJacobian                       True
pfset Solver.Nonlinear.DerivativeEpsilon                 1e-2

pfset Solver.Linear.KrylovDimension                      10

pfset Solver.Linear.Preconditioner                       MGSemi
pfset Solver.Linear.Preconditioner.MGSemi.MaxIter        1
pfset Solver.Linear.Preconditioner.MGSemi.MaxLevels      100

##################################################
# Output Silo Files
##################################################
pfset Solver.WriteSiloSubsurfData       True
pfset Solver.WriteSiloPressure          True
pfset Solver.WriteSiloSaturation        True
pfset Solver.WriteSiloConcentration     True
pfset Solver.WriteSiloMask              True  
pfset Solver.WriteSiloSlopes            True
pfset Solver.WriteSiloSpecificStorage   True
pfset Solver.WriteSiloMannings          True  

##################################################
# Run Program, Write Files
##################################################
# This part of the script converts the pressure pfb into a
# hydraulic head value

file mkdir "../parflow_Fbm"
cd  "../parflow_Fbm"
set inFile [pfload "../madepermFbm_slim/madepermFbm.sa"]
pfsave $inFile -pfb "madepermFbm.pfb"

pfdist madepermFbm.pfb
pfrun MADE
pfundist MADE
pfundist madepermFbm.pfb

# In order to have a pfb file for slim we are using indicator files for the vanGenuctin parameters
#set vgn [pfload "../vgn_homog.sa"] 
#pfsave $vgn -pfb "vgn_homog.pfb"
#set vga [pfload "../vga_homog.sa"] 
#pfsave $vga -pfb "vga_homog.pfb"
#set sres [pfload "../sres_homog.sa"] 
#pfsave $sres -pfb "sres_homog.pfb"
#set ssat [pfload "../ssat_homog.sa"] 
#pfsave $ssat -pfb "ssat_homog.pfb"


#set inFile [pfload "../satindicatorFile_lowKb.sa"]
#pfsave $inFile -pfb "satindicatorFile_lowKb.pfb"

#pfdist vgn_homog.pfb
#pfdist vga_homog.pfb
#pfdist sres_homog.pfb
#pfdist ssat_homog.pfb
#pfrun homog_MADE
#pfundist homog_MADE
#pfundist vgn_homog.pfb
#pfundist vga_homog.pfb
#pfundist sres_homog.pfb
#pfundist ssat_homog.pfb


#for {set k 0} {$k <= 460} {incr k 1} { 
#   set filename [format "homog_MADE.out.press.%05d.pfb" $k]
#   set press [pfload $filename]
#   set head [pfhhead $press] 
#   set filename [format "homog_MADE.out.head.%05d.silo" $k]
#  pfsave $head -silo $filename
#   }

#set press [pfload homog_MADE.out.press.pfb]
#set head [pfhhead $press]
#pfsave $head -silo homog.MADE.head.silo
#pfsave $head -pfb homog.MADE.head.pfb

cd ".."
