import ../measuremancer
import std / unittest
import std / math

proc `=~=`[T](a, b: Measurement[T]): bool =
  # echo a.value.float, " ± ", a.error.float, " and ", b.value.float, " ± ", b.error.float
  result = almostEqual(a.value.float, b.value.float) and almostEqual(a.error.float, b.error.float)

suite "Measurement & constant":
  test "Addition":
    let x = 1.0 ± 0.1
    let y = 5.0
    check x + y == 6.0 ± 0.1
    check y + x == 6.0 ± 0.1

  test "Subtraction":
    let x = 1.0 ± 0.1
    let y = 5.0
    check x - y == -4.0 ± 0.1
    check y - x == 4.0 ± 0.1

  test "Multiplication":
    let x = 1.0 ± 0.1
    let y = 5.0
    check x * y == 5.0 ± 0.5
    check y * x == 5.0 ± 0.5

  test "Division":
    let x = 5.0 ± 0.5
    let y = 5.0
    check x / y == 1.0 ± 0.1
    check y / x == 1.0 ± 0.1

  test "Division by int literal / const with Measurement":
    let x = 5.0 ± 0.5
    const y = 5
    check x / y == 1.0 ± 0.1
    check y / x == 1.0 ± 0.1

suite "Measurement & measurement":
  test "Addition":
    let x = 1.0 ± 0.1
    let y = 5.0 ± 0.5
    # This should be:
    # √( 0.1² + 0.5² )
    check x + y =~= 6.0 ± 0.5099019513592785
    check y + x =~= 6.0 ± 0.5099019513592785

  test "Subtraction":
    let x = 1.0 ± 0.1
    let y = 5.0 ± 0.5
    # This should be:
    # √( 0.1² + 0.5² )
    check x - y =~= -4.0 ± 0.5099019513592785
    check y - x =~=  4.0 ± 0.5099019513592785

  test "Multiplication":
    let x = 1.0 ± 0.1
    let y = 5.0 ± 0.5
    # This should be:
    # √( 1.0² · 0.5² + 5.0² · 0.1² ) = √ 0.5
    check x * y == 5.0 ± sqrt(0.5)
    check y * x == 5.0 ± sqrt(0.5)

  test "Division":
    let x = 1.0 ± 0.1
    let y = 5.0 ± 0.5
    # This should be:
    # √( 0.1² / 5.0² + 1.0² · 0.5² / 5.0⁴ ) = √ 8e-4 = 0.0282842712475
    check x / y =~= 0.2 ± 0.02828427124746191
    # This should be:
    # √( 0.5² / 1.0² + 5.0² · 0.1² / 1.0⁴ ) = √ 0.5 = 0.707106781187
    check y / x =~= 5.0 ± 0.7071067811865476

  test "Power (integer literal)":
    let x = 2.0 ± 0.2
    # This should be:
    # √( (2 · 2.0)² · 0.2² ) = 0.8
    check x ^ 2 == 4.0 ± 0.8
    check x ** 2 == 4.0 ± 0.8

  test "Power (static integer)":
    let x = 2.0 ± 0.2
    # This should be:
    # √( (2 · 2.0)² · 0.2² ) = 0.8
    const y = 2
    check x ^ y == 4.0 ± 0.8
    check x ** y == 4.0 ± 0.8

  test "Power (integer)":
    let x = 2.0 ± 0.2
    # This should be:
    # √( (2 · 2.0)² · 0.2² ) = 0.8
    let y = 2
    check x ^ y == 4.0 ± 0.8
    check x ** y == 4.0 ± 0.8

  test "Power (float)":
    let x = 2.0 ± 0.2
    # This should be:
    # √( [ 2.5 · 2.0^{1.5} ]² · 0.2² ) = √ ( 50 · 0.² ) = sqrt(2)
    check x ^ 2.5 == pow(2.0, 2.5) ± sqrt(2.0)
    check x ** 2.5 == pow(2.0, 2.5) ± sqrt(2.0)

  test "Negative power (integer)":
    let x = 2.0 ± 0.2
    # This should be:
    # √( (-2 / x³)² · Δx² ) = √( (-2 / 2³)² · 0.2² ) = 0.05
    check x ^ (-2) == 1.0 / (x^2)
    check x ** -2 == 1.0 / (x^2)
    check x ^ -2 == 0.250 ± 0.0500

  test "Power of two measurements":
    let x = 2.0 ± 0.2
    let y = 5.0 ± 0.5
    check x ^ y == 32.0 ± 19.467818870203708
    check x ** y == 32.0 ± 19.467818870203708

suite "Measurement and same symbol":
  # Simple correlations due to same symbol being used are correctly handled.
  # i.e. subtraction and division results in a measurement without error.
  test "Addition":
    let x = 1.0 ± 0.1
    check x + x =~= 2.0 ± 0.2

  test "Subtraction":
    let x = 1.0 ± 0.1
    check x - x == 0.0 ± 0.0

  test "Multiplication":
    let x = 1.0 ± 0.1
    check x * x =~= 1.0 ± 0.2

  test "Division":
    let x = 1.0 ± 0.1
    check x / x =~= 1.0 ± 0.0

suite "Trigonometric functions of measurements":
  test "sin":
    let x = 1.0 ± 0.1
    # This should be:
    # √ (cos²(1.0) · 0.1²)
    check sin(x) =~= sin(1.0) ± sqrt(cos(1.0)^2 * 0.1^2) # 0.05403023058681398

  test "cos":
    let x = 1.0 ± 0.1
    # This should be:
    # √ (sin²(1.0) · 0.1²)
    check cos(x) =~= cos(1.0) ± sqrt(sin(1.0)^2 * 0.1^2)

  test "tan":
    let x = 1.0 ± 0.1
    # This should be:
    # √ ( (1.0 + tan(1.0)²)² · 0.1² )
    check tan(x) =~= tan(1.0) ± sqrt(( (1.0 + tan(1.0)^2)^2 * 0.1^2) )

  ## And so on. The general error for these functions is simply the
  ## √ (∂f/∂x² Δx²)
  ## where the derivative can be looked up in the table in the main source file.
  test "sinh":
    let x = 1.0 ± 0.1
    check sinh(x) =~= sinh(1.0) ± sqrt(( cosh(1.0)^2 * 0.1^2) )

  test "cosh":
    let x = 1.0 ± 0.1
    check cosh(x) =~= cosh(1.0) ± sqrt(( sinh(1.0)^2 * 0.1^2) )

  test "arccos":
    let x = 0.5 ± 0.1
    check arccos(x) =~= arccos(0.5) ± sqrt( (-1.0 / sqrt(1.0 - (0.5^2)))^2 * 0.1^2 )

suite "More complicated examples":
  test "Function with measurement computing a gaussian":
    proc foo[T](x: T, μ, σ: float): T =
      let val = 2 * σ * σ
      result = exp(- (x*x) / val ) #/ (2 * σ^2))

    let x1 = 5.0 ± 1.0
    check foo(x1, 33.333, 100.0) =~= 0.9987507809245809 ± 0.0004993753904622905
    check foo(x1, 33.333, 100.0) =~= exp(- x1*x1 / (2 * 100.0 * 100.0))

import unchained
suite "Measurements of other types (e.g. unchained)":
  test "Addiiton of unchained units":
    let k1 = 5.0.keV ± 1.0.keV
    let k2 = 2.5.keV ± 1.5.keV
    check k1 + k2 =~= 7.5.keV ± 1.802775637731995.keV
    check( (k1 + k2).value.type is keV )
    check( (k1 + k2).error.type is keV )

  test "Construction of `measurement` via proc":
    let k1 = measurement(5.0.keV, 1.0.keV)
    let k2 = measurement(2.5.keV, 1.5.keV)
    check k1 + k2 =~= 7.5.keV ± 1.802775637731995.keV
    check( (k1 + k2).value.type is keV )
    check( (k1 + k2).error.type is keV )

  test "Addition of incompatible units fails":
    let k1 = 5.0.keV ± 1.0.keV
    let m = 1.0.kg ± 0.1.kg
    when compiles(k1 + m):
      error("This must not compile!")
    else:
      check true

  test "Multiplication yields new units":
    let m = 1.0.kg ± 0.1.kg
    let a = 9.81.m•s⁻² ± 0.5.m•s⁻²
    check m * a =~= 9.81.kg•m•s⁻² ± 1.101072658819571.kg•m•s⁻²
    check( (m * a).value.type is kg•m•s⁻² )
    check( (m * a).error.type is kg•m•s⁻² )

  # integer powers can be supported for units. Float however not! (floats representing
  # integer powers need to be converted to int manually!)
  test "Power (integer literal)":
    let x = 2.0.m ± 0.2.m
    # This should be:
    # √( (2 · 2.0)² · 0.2² ) = 0.8
    check x ^ 2 == 4.0.m² ± 0.8.m²
    check x ** 2 == 4.0.m² ± 0.8.m²

  test "Power (static integer)":
    let x = 2.0.m ± 0.2.m
    # This should be:
    # √( (2 · 2.0)² · 0.2² ) = 0.8
    const y = 2
    check x ^ y == 4.0.m² ± 0.8.m²
    check x ** y == 4.0.m² ± 0.8.m²

  test "Negative power with units":
    let x = 2.0.m ± 0.2.m
    defUnit(Meter⁻²)
    # This should be:
    # √( (-2 / x³)² · Δx² ) = √( (-2 / 2³)² · 0.2² ) = 0.05
    check x ^ -2 == 1.0 / (x^2)
    ## NOTE: As `**` does not use `unchained's` rewrite `^`, but rather
    ## the regular `**` proc, the resulting unit is equivalent, but not
    ## the *same Nim unit* (as they are `distinct` and `Meter⁻²` is not
    ## declared in the measuremancer context. As such we cannot check the
    ## sanity of the types at CT.
    # TODO: investigate why this prints as `0.0500` instead of `0.050`
    check $(x ** (-2)) == "0.250 ± 0.0500 Meter⁻²"
    check (x ** -2).to(Meter⁻²) == to(1.0 / (x^2), Meter⁻²)

  test "Trigonometric functions only work with UnitLess":
    block UnitLessCheck:
      ## the division is just a simple way to get a `UnitLess` value
      ## as on `unchained` HEAD (2022-11-16) `UnitLess` is not exported.
      let x = (1.0.kg / 1.0.kg) ± (0.1.kg / 1.0.kg)
      # This should be:
      # √ (cos²(1.0) · 0.1²)
      check typeof(x.value) is UnitLess
      # application of `sin` returns a `float`! Valid units for `sin` are those convertible
      # to `float` via converter defined in unchained
      check typeof(sin(x)) is Measurement[float]
      let exp = sin(x.value) ± sqrt(cos(x.value)^2 * (x.error)^2)
      check typeof(exp) is Measurement[float]
      check sin(x) =~= exp
    block RadianCheck: # converts to UnitLess
      let x = 1.0.Radian ± 0.1.Radian
      # This should be:
      # √ (cos²(1.0) · 0.1²)
      check sin(x) =~= sin(1.0.Radian) ± sqrt(cos(1.0.Radian)^2 * 0.1.Radian^2) # 0.05403023058681398
