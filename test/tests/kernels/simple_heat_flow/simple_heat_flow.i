[Mesh]
  [./my_mesh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 50
        ny = 5
        xmin = 0
        xmax = 0.1
        ymin = 0
        ymax = 0.01
    [../]
[]

[Variables]
  [./T]
        order = FIRST
        family = LAGRANGE
        initial_condition = 300 # K
  [../]
[]

[AuxVariables]
  # Parameters for Steel
  [./rho]
      order = FIRST
      family = LAGRANGE
      initial_condition = 7750  # kg/m^3
  [../]
  
  [./cp]
      order = FIRST
      family = LAGRANGE
      initial_condition = 466  # J/kg/K
  [../]
  
  [./K]
      order = FIRST
      family = LAGRANGE
      initial_condition = 45   # W/m/K
  [../]
  
  [./S]
      order = FIRST
      family = LAGRANGE
      initial_condition = 1e6   # W/m^3
  [../]
  
  [./Tamb]
      order = FIRST
      family = LAGRANGE
      initial_condition = 273   # W/m^3
  [../]
  
[]

[Kernels]
  [./heat_accum]
    type = HeatAccumulation
    variable = T
	density = rho
	heat_capacity = cp
  [../]
  [./heat_cond]
    type = HeatConduction
    variable = T
	thermal_conductivity = K
  [../]
  [./heat_conv]
    type = HeatConvection
    variable = T
	coupled_temperature = Tamb
	convection_coeff = 2e2
	specific_area = 2e2
  [../]
  [./heat_source]
    type = HeatSource
    variable = T
	coupled_source = S
  [../]
[]

[BCs]
[./left]
    type = DirichletBC
    variable = T
    boundary = 'left'
    value = 350
  [../]

[]

[Postprocessors]	

	[./T_left]
        type = SideAverageValue
        boundary = 'left'
        variable = T
        execute_on = 'initial timestep_end'
    [../]
 
    [./T_right]
        type = SideAverageValue
        boundary = 'right'
        variable = T
        execute_on = 'initial timestep_end'
    [../]
	
	[./T_avg]
      type = ElementAverageValue
      # block = NAME_OF_SUBDOMAIN  # Optional if block has different names
      variable = T
      execute_on = 'initial timestep_end'
  [../]
[]

[Preconditioning]
    [./SMP_PJFNK]
      type = SMP
      full = true
      solve_type = pjfnk
    [../]

[] #END Preconditioning

[Executioner]
  type = Transient
  scheme = implicit-euler
  
  start_time = 0.0
  end_time = 1000.0
  dtmax = 100.0

  [./TimeStepper]
    type = ConstantDT
    dt = 100.0
  [../]
  
  petsc_options = '-snes_converged_reason

                      -ksp_gmres_modifiedgramschmidt'

    # NOTE: The sub_pc_type arg not used if pc_type is ksp,
    #       Instead, set the ksp_ksp_type to the pc method
    #       you want. Then, also set the ksp_pc_type to be
    #       the terminal preconditioner.
    #
    # Good terminal precon options: lu, ilu, asm, gasm, pbjacobi
    #                               bjacobi, redundant, telescope
    petsc_options_iname ='-ksp_type
                          -pc_type

                          -sub_pc_type

                          -snes_max_it

                          -sub_pc_factor_shift_type
                          -pc_asm_overlap

                          -snes_atol
                          -snes_rtol

                          -ksp_ksp_type
                          -ksp_pc_type'

    # snes_max_it = maximum non-linear steps
    petsc_options_value = 'fgmres
                           ksp

                           ilu

                           10
                           NONZERO
                           10
                           1E-8
                           1E-10

                           gmres
                           ilu'

    line_search = none   # none, bt, l2, basic
    nl_rel_tol = 1e-8
    nl_abs_tol = 1e-8
    nl_rel_step_tol = 1e-12
    nl_abs_step_tol = 1e-12
    nl_max_its = 10
    l_tol = 1e-6
    l_max_its = 300
[]

[Outputs]
  print_linear_residuals = true
  exodus = true
  csv = true
[]
