program ICE_docgen
    use ICE_Read_Sim_Param, only: Register_Sim_Param
    use ICE_Read_IO,        only: Register_IO_Fields, Register_Probes
    use ICE_Read_Numerics,  only: Register_Numerics, Register_Families
    use ICE_Read_Physics,   only: Register_Physics
    use ICE_Input_Registry
    implicit none

    ! Build registry entries
    call Register_Sim_Param()
    call Register_IO_Fields()
    call Register_Probes(1, 'probe-section')
    call Register_Numerics(2)
    call Register_Physics()
    call Register_Families(1)

    call reg%generate_markdown('docs/user/registry.md')

end program ICE_docgen
