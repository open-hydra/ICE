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
    logical :: unread       = .false.   ! in input.ini with a value that does not read as its type (Load_Ini)
    character(len=:), allocatable :: text   ! the value as written in input.ini, when the key is there (Load_Ini)
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

  !> The first error of the values read, '' when there is none. A key whose value does not read
  !  as its type (Load_Ini) comes first, then a required key not set, then a real given in
  !  input.ini that is not finite (every comparison with a NaN is false and an infinity passes
  !  a bound such as > 0; gfortran reads 1e999 as an infinity), then the rule of the column
  !  Allowed: a bound (> x, >= x, < x, <= x) or an interval (in (a, b), [ or ] for an end that
  !  belongs to it) on a number and on each value of a list of numbers, a list of names for a
  !  string, logical for a logical. A rule of another form is an error of the registry itself,
  !  so that no rule is skipped in silence.
  function Validate_Registry() result(out)
    use ICE_Finite_m, only: Not_Finite
    character(len=1024) :: out
    integer  :: i, k
    real(R8) :: val
    character(len=16) :: nth
    out = ""
    do i = 1, reg%size
      if (reg%params(i)%unread) then
        out = Unread_Message(i)
        return
      end if
      if (reg%params(i)%required .and. .not. reg%params(i)%is_set) then
        out = "[ERROR] Required parameter not set: "//trim(reg%params(i)%name)
        return
      end if
      ! An unset key keeps its default; the values of a per-material list are Setup_Materials'
      if (reg%params(i)%is_set .and. .not. reg%params(i)%multi) then
        if (associated(reg%params(i)%value%r)) then
          if (Not_Finite(reg%params(i)%value%r)) then
            out = "[ERROR] "//Key_And_Value(i)//": the value is not a finite real number"
            return
          end if
        else if (associated(reg%params(i)%value%rarr)) then
          do k = 1, size(reg%params(i)%value%rarr)
            if (Not_Finite(reg%params(i)%value%rarr(k))) then
              write(nth,'(I0)') k
              out = "[ERROR] "//Key_And_Value(i)//": value "//trim(nth)//" of the list is not a finite real number"
              return
            end if
          end do
        end if
      end if
      if (reg%params(i)%allowed == "") cycle
      ! A list of numbers: the rule holds for each of its values
      if (associated(reg%params(i)%value%iarr)) then
        do k = 1, size(reg%params(i)%value%iarr)
          val = real(reg%params(i)%value%iarr(k), R8)
          call validate_numeric(reg%params(i)%name, val, reg%params(i)%allowed, out)
          if (out /= "") exit
        end do
      else if (associated(reg%params(i)%value%rarr)) then
        do k = 1, size(reg%params(i)%value%rarr)
          val = reg%params(i)%value%rarr(k)
          call validate_numeric(reg%params(i)%name, val, reg%params(i)%allowed, out)
          if (out /= "") exit
        end do
      end if
      if (associated(reg%params(i)%value%iarr) .or. associated(reg%params(i)%value%rarr)) then
        if (out /= "") then
          write(nth,'(I0)') k
          out = trim(out)//" (value "//trim(nth)//" of the list)"
          return
        end if
        cycle
      end if
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
      case(TYPE_LOG)
        ! FiNeR takes a logical only as a logical (Load_Ini marks any other value): nothing is
        ! left to check, if the rule says no more than that
        if (.not. logical_rule(reg%params(i)%allowed)) then
          out = "[ERROR] "//trim(reg%params(i)%name)//": the registry does not check its rule "// &
                trim(reg%params(i)%allowed)
          return
        end if
      case default
        cycle
      end select
    end do
  end function Validate_Registry

  !> The key i of input.ini as the messages give it: the file, the section, the key and the
  !  value as written, when the file holds one.
  function Key_And_Value(i) result(out)
    integer, intent(in) :: i
    character(len=:), allocatable :: out
    out = "input.ini, ["//trim(reg%params(i)%section)//"] "//trim(reg%params(i)%name)//" ="
    if (allocated(reg%params(i)%text)) out = out//" "//trim(reg%params(i)%text)
  end function Key_And_Value

  !> The error of a key whose value does not read as the type of its parameter (Load_Ini).
  function Unread_Message(i) result(out)
    integer, intent(in) :: i
    character(len=:), allocatable :: out
    character(len=16) :: n
    out = "[ERROR] "//Key_And_Value(i)//": the value "
    if (.not. allocated(reg%params(i)%text)) then
      out = out//"is empty"
    else if (len_trim(reg%params(i)%text) == 0) then
      out = out//"is empty"
    else if (associated(reg%params(i)%value%iarr)) then
      write(n,'(I0)') size(reg%params(i)%value%iarr)
      out = out//"is not a list of at most "//trim(n)//" integers"
    else if (associated(reg%params(i)%value%rarr)) then
      write(n,'(I0)') size(reg%params(i)%value%rarr)
      out = out//"is not a list of at most "//trim(n)//" finite real numbers"
    else
      select case (reg%params(i)%type_id)
      case (TYPE_INT)
        out = out//"is not an integer"
      case (TYPE_REAL)
        out = out//"is not a finite real number"
      case (TYPE_LOG)
        out = out//"is not a logical (true or false)"
      case default
        out = out//"cannot be read"
      end select
    end if
  end function Unread_Message

  !> The rule of a logical: logical, or the list of its two values (with or without the periods).
  logical function logical_rule(rule)
    character(*), intent(in) :: rule
    character(len=:), allocatable :: r
    integer :: k
    r = ""
    do k = 1, len(rule)
      if (rule(k:k) /= " " .and. rule(k:k) /= ".") r = r//rule(k:k)
    end do
    logical_rule = r == "logical" .or. r == "true,false" .or. r == "false,true"
  end function logical_rule

  !> The rule of a number, read from its start: a bound (> x, >= x, < x, <= x) or an interval,
  !  in (a, b), where ( and ) leave the end out and [ and ] take it in. A rule of another form,
  !  or whose numbers do not read, is an error of the registry itself.
  subroutine validate_numeric(name, val, rule, out)
    character(*), intent(in)    :: name, rule
    real(R8),     intent(in)    :: val
    character(*), intent(inout) :: out
    real(R8) :: limit, lo, hi
    character(len=:), allocatable :: r
    integer :: c, n, ios1, ios2
    logical :: ok, known
    r = trim(adjustl(rule))
    known = .true.
    ok    = .true.
    ios1  = 0
    ios2  = 0
    if (index(r,"in") == 1) then
      ! the first and last characters after "in" open and close the interval
      r = trim(adjustl(r(3:)))
      n = len(r)
      c = index(r,",")
      known = n >= 5 .and. c > 2 .and. c < n - 1
      if (known) known = scan(r(1:1),"([") == 1 .and. scan(r(n:n),")]") == 1
      if (known) then
        read(r(2:c-1),*,iostat=ios1) lo
        read(r(c+1:n-1),*,iostat=ios2) hi
        known = ios1 == 0 .and. ios2 == 0
      end if
      if (known) ok = merge(val >= lo, val > lo, r(1:1) == "[") .and. merge(val <= hi, val < hi, r(n:n) == "]")
      r = "in "//r
    else if (index(r,">=") == 1) then
      read(r(3:),*,iostat=ios1) limit
      ok = val >= limit
    else if (index(r,"<=") == 1) then
      read(r(3:),*,iostat=ios1) limit
      ok = val <= limit
    else if (index(r,">") == 1) then
      read(r(2:),*,iostat=ios1) limit
      ok = val > limit
    else if (index(r,"<") == 1) then
      read(r(2:),*,iostat=ios1) limit
      ok = val < limit
    else
      known = .false.
    end if
    if (.not. known .or. ios1 /= 0) then
      out = "[ERROR] "//trim(name)//": the registry does not check its rule "//trim(rule)
    else if (.not. ok) then
      out = "[ERROR] "//trim(name)//" must be "//r
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
