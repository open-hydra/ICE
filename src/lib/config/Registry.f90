module ICE_Input_Registry

  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private

  integer, parameter :: TYPE_INT=1, TYPE_REAL=2, TYPE_LOG=3, TYPE_STR=4

  public :: registry_t, Validate_Registry, param_is_set

  !--------------------------------------------------------
  ! Value container (typed pointers)
  !--------------------------------------------------------
  type :: param_value_t
    integer,          pointer :: i    => null()
    real(R8),         pointer :: r    => null()
    logical,          pointer :: l    => null()
    character(len=:), pointer :: s    => null()
    integer,          pointer :: iarr(:) => null()
    real(R8),         pointer :: rarr(:) => null()
  end type

  !--------------------------------------------------------
  ! Parameter description
  !--------------------------------------------------------
  type :: param_t
    character(len=:), allocatable :: section
    character(len=:), allocatable :: name
    character(len=:), allocatable :: description
    character(len=:), allocatable :: default_str
    character(len=:), allocatable :: allowed
    logical :: required = .false.
    logical :: is_set   = .false.
    logical :: per_material = .false.   ! may hold one value per material (Setup_Materials reads them)
    logical :: multi        = .false.   ! holds several values: the scalar target keeps its default
    integer :: type_id  = 0
    type(param_value_t) :: value
  end type

  !--------------------------------------------------------
  ! Registry container
  !--------------------------------------------------------
  type :: registry_t
    type(param_t), allocatable :: params(:)
    integer :: size     = 0
    integer :: capacity = 0
  contains
    procedure :: reserve
    procedure :: add_int
    procedure :: add_real
    procedure :: add_logical
    procedure :: add_string
    procedure :: add_int_array
    procedure :: add_real_array
    generic   :: add => add_int, add_real, add_logical, add_string, &
                        add_int_array, add_real_array
    procedure :: generate_markdown
  end type

  type(registry_t), public :: reg

contains

  subroutine reserve(this, newcap)
    class(registry_t), intent(inout) :: this
    integer,           intent(in)    :: newcap
    type(param_t), allocatable :: tmp(:)
    if (newcap <= this%capacity) return
    allocate(tmp(newcap))
    if (this%size > 0) tmp(1:this%size) = this%params(1:this%size)
    call move_alloc(tmp, this%params)
    this%capacity = newcap
  end subroutine reserve

  subroutine ensure_space(this)
    class(registry_t), intent(inout) :: this
    if (this%size == this%capacity) then
      if (this%capacity == 0) then
        call this%reserve(8)
      else
        call this%reserve(2*this%capacity)
      end if
    end if
  end subroutine ensure_space

  subroutine add_int(this, section, name, var, default, desc, allowed, required)
    class(registry_t), intent(inout) :: this
    character(*),      intent(in)    :: section, name, default, desc, allowed
    logical,           intent(in)    :: required
    integer, target,   intent(inout) :: var
    integer :: n
    call ensure_space(this)
    this%size = this%size + 1 ; n = this%size
    this%params(n)%section     = section
    this%params(n)%name        = name
    this%params(n)%description = desc
    this%params(n)%default_str = default
    this%params(n)%allowed     = allowed
    this%params(n)%required    = required
    this%params(n)%type_id     = TYPE_INT
    this%params(n)%value%i    => var
    read(default,*) var
  end subroutine add_int

  subroutine add_real(this, section, name, var, default, desc, allowed, required, per_material)
    class(registry_t), intent(inout) :: this
    character(*),      intent(in)    :: section, name, default, desc, allowed
    logical,           intent(in)    :: required
    real(R8), target,  intent(inout) :: var
    logical, optional, intent(in)    :: per_material
    integer :: n
    call ensure_space(this)
    this%size = this%size + 1 ; n = this%size
    this%params(n)%section     = section
    this%params(n)%name        = name
    this%params(n)%description = desc
    this%params(n)%default_str = default
    this%params(n)%allowed     = allowed
    this%params(n)%required    = required
    this%params(n)%type_id     = TYPE_REAL
    this%params(n)%value%r    => var
    if (present(per_material)) this%params(n)%per_material = per_material
    read(default,*) var
  end subroutine add_real

  subroutine add_logical(this, section, name, var, default, desc, allowed, required)
    class(registry_t), intent(inout) :: this
    character(*),      intent(in)    :: section, name, default, desc, allowed
    logical,           intent(in)    :: required
    logical, target,   intent(inout) :: var
    integer :: n
    call ensure_space(this)
    this%size = this%size + 1 ; n = this%size
    this%params(n)%section     = section
    this%params(n)%name        = name
    this%params(n)%description = desc
    this%params(n)%default_str = default
    this%params(n)%allowed     = allowed
    this%params(n)%required    = required
    this%params(n)%type_id     = TYPE_LOG
    this%params(n)%value%l    => var
    read(default,*) var
  end subroutine add_logical

  subroutine add_string(this, section, name, var, default, desc, allowed, required)
    class(registry_t), intent(inout) :: this
    character(*),      intent(in)    :: section, name, default, desc, allowed
    logical,           intent(in)    :: required
    character(len=*), target, intent(inout) :: var
    integer :: n
    call ensure_space(this)
    this%size = this%size + 1 ; n = this%size
    this%params(n)%section     = section
    this%params(n)%name        = name
    this%params(n)%description = desc
    this%params(n)%default_str = default
    this%params(n)%allowed     = allowed
    this%params(n)%required    = required
    this%params(n)%type_id     = TYPE_STR
    this%params(n)%value%s    => var
    var = default
  end subroutine add_string

  subroutine add_int_array(this, section, name, var, default, desc, allowed, required)
    class(registry_t), intent(inout) :: this
    character(*),      intent(in)    :: section, name, default, desc, allowed
    logical,           intent(in)    :: required
    integer, target,   intent(inout) :: var(:)
    integer :: n, defval
    call ensure_space(this)
    this%size = this%size + 1 ; n = this%size
    this%params(n)%section        = section
    this%params(n)%name           = name
    this%params(n)%description    = desc
    this%params(n)%default_str    = default
    this%params(n)%allowed        = allowed
    this%params(n)%required       = required
    this%params(n)%type_id        = TYPE_INT
    this%params(n)%value%iarr    => var
    read(default,*) defval ; var(:) = defval
  end subroutine add_int_array

  subroutine add_real_array(this, section, name, var, default, desc, allowed, required)
    class(registry_t), intent(inout) :: this
    character(*),      intent(in)    :: section, name, default, desc, allowed
    logical,           intent(in)    :: required
    real(R8), target,  intent(inout) :: var(:)
    integer  :: n
    real(R8) :: defval
    call ensure_space(this)
    this%size = this%size + 1 ; n = this%size
    this%params(n)%section        = section
    this%params(n)%name           = name
    this%params(n)%description    = desc
    this%params(n)%default_str    = default
    this%params(n)%allowed        = allowed
    this%params(n)%required       = required
    this%params(n)%type_id        = TYPE_REAL
    this%params(n)%value%rarr    => var
    read(default,*) defval ; var(:) = defval
  end subroutine add_real_array

  !> True when the INI gave a value for this key (false for a key never registered).
  logical function param_is_set(section, name)
    character(*), intent(in) :: section, name
    integer :: i
    param_is_set = .false.
    do i = 1, reg%size
      if (reg%params(i)%section == section .and. reg%params(i)%name == name) then
        param_is_set = reg%params(i)%is_set
        return
      end if
    end do
  end function param_is_set

  !> True when another entry with the same real target (a key and its alias) was set by the INI.
  logical function set_by_alias(i)
    integer, intent(in) :: i
    integer :: j
    set_by_alias = .false.
    do j = 1, reg%size
      if (j == i .or. .not. reg%params(j)%is_set) cycle
      if (.not. associated(reg%params(j)%value%r)) cycle
      if (associated(reg%params(j)%value%r, reg%params(i)%value%r)) then
        set_by_alias = .true.
        return
      end if
    end do
  end function set_by_alias

  function Validate_Registry() result(out)
    character(len=1024) :: out
    integer  :: i
    real(R8) :: val
    out = ""
    do i = 1, reg%size
      if (reg%params(i)%required .and. .not. reg%params(i)%is_set) then
        out = "[ERROR] Required parameter not set: "//trim(reg%params(i)%name)
        return
      end if
      if (reg%params(i)%allowed == "") cycle
      if (associated(reg%params(i)%value%iarr) .or. associated(reg%params(i)%value%rarr)) cycle
      select case(reg%params(i)%type_id)
      case(TYPE_INT)
        if (.not. associated(reg%params(i)%value%i)) cycle
        val = real(reg%params(i)%value%i, R8)
        call validate_numeric(reg%params(i)%name, val, reg%params(i)%allowed, out)
        if (out /= "") return
      case(TYPE_REAL)
        if (.not. associated(reg%params(i)%value%r)) cycle
        if (reg%params(i)%multi) cycle   ! each value is checked by Setup_Materials
        if (set_by_alias(i)) cycle   ! checked under the name given; Setup_Materials refuses both names
        val = reg%params(i)%value%r
        call validate_numeric(reg%params(i)%name, val, reg%params(i)%allowed, out)
        if (out /= "") return
      case(TYPE_STR)
        if (.not. associated(reg%params(i)%value%s)) cycle
        if (trim(reg%params(i)%value%s) == "") cycle
        call validate_string(reg%params(i)%name, reg%params(i)%value%s, reg%params(i)%allowed, out)
        if (out /= "") return
      case default
        cycle
      end select
    end do
  end function Validate_Registry

  subroutine validate_numeric(name, val, rule, out)
    character(*), intent(in)    :: name, rule
    real(R8),     intent(in)    :: val
    character(*), intent(inout) :: out
    real(R8) :: limit
    if (index(rule,">=") > 0) then
      read(rule(index(rule,">=")+2:),*) limit
      if (val < limit) out = "[ERROR] "//trim(name)//" must be "//trim(rule)
    else if (index(rule,"<=") > 0) then
      read(rule(index(rule,"<=")+2:),*) limit
      if (val > limit) out = "[ERROR] "//trim(name)//" must be "//trim(rule)
    else if (index(rule,">") > 0) then
      read(rule(index(rule,">")+1:),*) limit
      if (val <= limit) out = "[ERROR] "//trim(name)//" must be "//trim(rule)
    else if (index(rule,"<") > 0) then
      read(rule(index(rule,"<")+1:),*) limit
      if (val >= limit) out = "[ERROR] "//trim(name)//" must be "//trim(rule)
    end if
  end subroutine validate_numeric

  subroutine validate_string(name, value, allowed, out)
    character(*), intent(in)    :: name, value, allowed
    character(*), intent(inout) :: out
    character(len=:), allocatable :: token, lower_value
    integer :: i, istart, iend
    logical :: is_sep
    istart = 1
    do i = 1, len(allowed)+1
      is_sep = i > len(allowed)
      if (.not. is_sep) is_sep = (allowed(i:i) == ",")
      if (is_sep) then
        iend  = i - 1
        token = adjustl(allowed(istart:iend))
        lower_value = value
        if (trim(token) == trim(lower_value)) return
        istart = i + 1
      end if
    end do
    out = "[ERROR] "//trim(name)//" must be one of: "//trim(allowed)
  end subroutine validate_string

  subroutine generate_markdown(this, filename)
    class(registry_t), intent(in)           :: this
    character(*),      intent(in), optional :: filename
    integer :: i, j, unit
    character(len=:), allocatable :: fileout, current_section
    logical :: new_section
    if (present(filename)) then
      fileout = filename
    else
      fileout = "registry.md"
    end if
    open(newunit=unit, file=fileout, status="replace")
    write(unit,'(A)') "# Input Parameters"
    write(unit,'(A)') ""
    write(unit,'(A)') "Generated from the input registry by `bin/DocGen`; regenerate it after"
    write(unit,'(A)') "changing any `reg%add` call. Every parameter is optional unless the"
    write(unit,'(A)') "Required column says otherwise, and omitting one selects the default."
    !> One table per section, gathering entries wherever they were registered: the
    !  sections are filled by several Register_* routines, so consecutive entries do
    !  not all belong to the same one.
    do i = 1, this%size
      new_section = .true.
      do j = 1, i - 1
        if (trim(this%params(j)%section) == trim(this%params(i)%section)) then
          new_section = .false.
          exit
        end if
      end do
      if (.not. new_section) cycle
      current_section = this%params(i)%section
      write(unit,'(A)') ""
      write(unit,'(A)') "## ["//trim(current_section)//"]"
      write(unit,'(A)') ""
      write(unit,'(A)') "| Parameter | Default | Allowed | Required | Description |"
      write(unit,'(A)') "|-----------|---------|---------|----------|-------------|"
      do j = i, this%size
        if (trim(this%params(j)%section) /= trim(current_section)) cycle
        write(unit,'(A)') "| `"//trim(this%params(j)%name)// &
            "` | "//trim(this%params(j)%default_str)// &
            " | "//trim(this%params(j)%allowed)// &
            " | "//merge("yes"," no", this%params(j)%required)// &
            " | "//trim(this%params(j)%description)//" |"
      end do
    end do
    close(unit)
  end subroutine generate_markdown

end module ICE_Input_Registry
