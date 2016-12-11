# solves 1-D convection diffusion to produce figure showing 
# how solution becomes unstable for Pe > 1

[GlobalParams]
  vel_x = vel_x
  porosity = 1.0
[]

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
  xmin = 0.0
  xmax = 1.0
[]

[Variables]
  [./rho]
  [../]
[]

[AuxVariables]
  [./vel_x]
  [../]
[]

[Functions]
  [./vel_x]
    type = ParsedFunction
    value = 10
  [../]
[]

[AuxKernels]
  [./vel_x]
    type = FunctionAux
    function = 'vel_x'
    variable = vel_x
  [../]
[]

[Kernels]
  #active = 'continuity diffusion' # no stabilization
  [./continuity]
    type = ContinuityEqn
    variable = rho
  [../]
  [./continuitySUPG]
    type = ContinuityEqnSUPG
    variable = rho
  [../]
  [./diffusion]
    type = Diffusion
    variable = rho
  [../]
  [./diffusionSUPG] # zero contribution if rho is linear or lower order
    type = FluidEnergyDiffusionSUPG
    variable = rho
    k_fluid = 1.0
    T_fluid = rho
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = rho
    boundary = 'left'
    value = 0.0
  [../]
  [./right]
    type = DirichletBC
    variable = rho
    boundary = 'right'
    value = 1.0
  [../]
[]

[Materials]
  [./tau]
    type = Tau1DConvectionDiffusion
    k = 1
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
    exodus = true
  [./console]
    type = Console
    max_rows = 400
  [../]
[]
