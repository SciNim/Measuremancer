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

  test "Power":
    let x = 5.0 ± 0.5
    let y = 2.0
    check x ** y == 25.0 ± 5.0
    when compiles((block:
      check y ** x # TODO check expected result
    )): check true
    else: check false

  test "Log":
    let x = 5.0 ± 0.5
    let y = 2.0
    # something weird here
    # log(x) complains missing second argument (base)
    # log(x, y) complains extra argument given (position 2)
    when compiles((block:
      check log(x, y) # TODO check expected result
      check log(y, x) # TODO check expected result
    )): check true
    else: check false

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

  test "Power":
    let x = 5.0 ± 0.5
    let y = 2.0 ± 0.5
    when compiles((block:
      check x ** y # TODO check expected result
      check y ** x # TODO check expected result
    )): check true
    else: check false

  test "Log":
    let x = 5.0 ± 0.5
    let y = 2.0 ± 0.5
    # something weird here
    # log(x) complains missing second argument (base)
    # log(x, y) complains extra argument given (position 2)
    when compiles((block:
      check log(x, y) # TODO check expected result
      check log(y, x) # TODO check expected result
    )): check true
    else: check false

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

  test "Power":
    let x = 5.0 ± 0.5
    when compiles((block:
      check x ** x == 25.0 ± 5.0 # check expected result
    )): check true
    else: check false

  test "Log":
    let x = 5.0 ± 0.5
    # something weird here
    # log(x) complains missing second argument (base)
    # log(x, y) complains extra argument given (position 2)
    when compiles((block:
      check log(x, x) # TODO check expected result
    )): check true
    else: check false

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

# can't access .val attribute, must rely on .value func
# this should be fixed in lib
proc `==`[T: FloatLike](m1, m2: Measurement[T]): bool =
  result = almostEqual(m1.value.float, m2.value.float) and almostEqual(m1.error.float, m2.error.float)

suite "[Unchained] Measurement & constant":
  test "Addition":
    let x = 1.0.m ± 0.1.m
    let y = 5.0.m
    check x + y == 6.0.m ± 0.1.m
    check y + x == 6.0.m ± 0.1.m

  test "Subtraction":
    let x = 1.0.m ± 0.1.m
    let y = 5.0.m
    check x - y == -4.0.m ± 0.1.m
    check y - x == 4.0.m ± 0.1.m

  test "Multiplication":
    let x = 1.0.m ± 0.1.m
    let y = 5.0.m
    check x * y == 5.0.m² ± 0.5.m²
    check y * x == 5.0.m² ± 0.5.m²

  test "Division":
    let x = 5.0.m ± 0.5.m
    let y = 5.0.m
    check x / y == 1.0.UnitLess ± 0.1.UnitLess
    check y / x == 1.0.UnitLess ± 0.1.UnitLess

  test "Power":
    let x = 5.0.m ± 0.5.m
    let y = 2.0.m
    when compiles((block:
      check x ** y == 25.0.m ± 5.0.m
      check y ** x == 1.0.m ± 0.1.m # check expected result
    )): check true
    else: check false

  test "Log":
    let x = 5.0.m ± 0.5.m
    let y = 2.0.m
    # something weird here
    # log(x) complains missing second argument (base)
    # log(x, y) complains extra argument given (position 2)
    when compiles((block:
      check log(x, y) # TODO check expected result
      check log(y, x) # TODO check expected result
    )): check true
    else: check false

suite "[Unchained] Measurement & measurement":
  test "Addition":
    let x = 1.0.m ± 0.1.m
    let y = 5.0.m ± 0.5.m
    # This should be:
    # √( 0.1² + 0.5² )
    check x + y =~= 6.0.m ± 0.5099019513592785.m
    check y + x =~= 6.0.m ± 0.5099019513592785.m

  test "Subtraction":
    let x = 1.0.m ± 0.1.m
    let y = 5.0.m ± 0.5.m
    # This should be:
    # √( 0.1² + 0.5² )
    check x - y =~= -4.0.m ± 0.5099019513592785.m
    check y - x =~=  4.0.m ± 0.5099019513592785.m

  test "Multiplication":
    let x = 1.0.m ± 0.1.m
    let y = 5.0.m ± 0.5.m
    # This should be:
    # √( 1.0² · 0.5² + 5.0² · 0.1² ) = √ 0.5
    check x * y == 5.0.m² ± sqrt(0.5).m²
    check y * x == 5.0.m² ± sqrt(0.5).m²

  test "Division":
    let x = 1.0.m ± 0.1.m
    let y = 5.0.m ± 0.5.m
    # This should be:
    # √( 0.1² / 5.0² + 1.0² · 0.5² / 5.0⁴ ) = √ 8e-4 = 0.0282842712475
    check x / y =~= 0.2.UnitLess ± 0.02828427124746191.UnitLess
    # This should be:
    # √( 0.5² / 1.0² + 5.0² · 0.1² / 1.0⁴ ) = √ 0.5 = 0.707106781187
    check y / x =~= 5.0.UnitLess ± 0.7071067811865476.UnitLess

  test "Power":
    let x = 5.0.m ± 0.5.m
    let y = 2.0.m ± 0.5.m
    when compiles((block:
      check x ** y == 25.0.m ± 5.0.m # check expected result
      check y ** x == 1.0.m ± 0.1.m # check expected result
    )): check true
    else: check false

  test "Log":
    let x = 5.0.m ± 0.5.m
    let y = 2.0.m ± 0.5.m
    # something weird here
    # log(x) complains missing second argument (base)
    # log(x, y) complains extra argument given (position 2)
    when compiles((block:
      check log(x, y) # TODO check expected result
      check log(y, x) # TODO check expected result
    )): check true
    else: check false

suite "[Unchained] Measurement and same symbol":
  # Simple correlations due to same symbol being used are correctly handled.
  # i.e. subtraction and division results in a measurement without error.
  test "Addition":
    let x = 1.0.m ± 0.1.m
    check x + x =~= 2.0.m ± 0.2.m

  test "Subtraction":
    let x = 1.0.m ± 0.1.m
    check x - x == 0.0.m ± 0.0.m

  test "Multiplication":
    let x = 1.0.m ± 0.1.m
    check x * x =~= 1.0.m² ± 0.2.m²

  test "Division":
    let x = 1.0.m ± 0.1.m
    check x / x =~= 1.0.UnitLess ± 0.0.UnitLess

  test "Power":
    let x = 5.0.m ± 0.5.m
    when compiles((block:
      check x ** x == 25.0.m ± 5.0.m # check expected result
    )): check true
    else: check false

  test "Log":
    let x = 5.0.m ± 0.5.m
    # something weird here
    # log(x) complains missing second argument (base)
    # log(x, y) complains extra argument given (position 2)
    when compiles((block:
      check log(x, x) # TODO check expected result
    )): check true
    else: check false

suite "[Unchained] Trigonometric functions of measurements":
  test "sin":
    let x = 1.0.m ± 0.1.m
    # This should be:
    # √ (cos²(1.0) · 0.1²)
    when compiles((block:
      check sin(x) =~= sin(1.0).m ± sqrt(cos(1.0)^2 * 0.1^2).m # 0.05403023058681398
    )): check true
    else: check false

  test "cos":
    let x = 1.0.m ± 0.1.m
    # This should be:
    # √ (sin²(1.0) · 0.1²)
    when compiles((block:
      check cos(x) =~= cos(1.0).m ± sqrt(sin(1.0)^2 * 0.1^2).m
    )): check true
    else: check false

  test "tan":
    let x = 1.0.m ± 0.1.m
    # This should be:
    # √ ( (1.0 + tan(1.0)²)² · 0.1² )
    when compiles((block:
      check tan(x) =~= tan(1.0).m ± sqrt(( (1.0 + tan(1.0)^2)^2 * 0.1^2) ).m
    )): check true
    else: check false

  ## And so on. The general error for these functions is simply the
  ## √ (∂f/∂x² Δx²)
  ## where the derivative can be looked up in the table in the main source file.
  test "sinh":
    let x = 1.0.m ± 0.1.m
    when compiles((block:
      check sinh(x) =~= sinh(1.0).m ± sqrt(( cosh(1.0)^2 * 0.1^2) ).m
    )): check true
    else: check false

  test "cosh":
    let x = 1.0.m ± 0.1.m
    when compiles((block:
      check cosh(x) =~= cosh(1.0).m ± sqrt(( sinh(1.0)^2 * 0.1^2) ).m
    )): check true
    else: check false

  test "arccos":
    let x = 0.5.m ± 0.1.m
    when compiles((block:
       check arccos(x) =~= arccos(0.5).m ± sqrt( (-1.0 / sqrt(1.0 - (0.5^2)))^2 * 0.1^2 ).m
    )): check true
    else: check false
  
  test "Function with measurement computing a gaussian":
    proc foo[T](x: T, μ, σ: float): T =
      let val = 2 * σ * σ
      result = exp(- (x*x) / val ) #/ (2 * σ^2))

    let x1 = 5.0.m ± 1.0.m
    when compiles((block:
      check foo(x1, 33.333, 100.0) =~= 0.9987507809245809 ± 0.0004993753904622905
      check foo(x1, 33.333, 100.0) =~= exp(- x1*x1 / (2 * 100.0 * 100.0))
    )): check true
    else: check false

  test "Addition of unchained units":
    let k1 = 5.0.keV ± 1.0.keV
    let k2 = 2.5.keV ± 1.5.keV
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
