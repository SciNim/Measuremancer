when not declared(assert): import std/assertions
import math, tables, strformat
from strutils import formatBiggestFloat, FloatFormatMode

export math
#[

This library follows the ideas of the excellent Measurements.jl package:
https://github.com/JuliaPhysics/Measurements.jl
and parts are essentially a straight up port of it (the types & idea of the `result`
procedure to compute the gradients in an elegant way)


This module ``must`` be compiled with `--experimental:unicodeOperators`!

]#

type
  FloatLike* = concept x
    x.toFloat is float

  ## A concept for any float like that is *not* an integer type
  FloatLikeNotInt* = concept x, type T
    x.toFloat is float
    not (T is SomeInteger)

  ## A concept for types that are float like and whose application of
  ## `pow` compiles and yields the same type. E.g. `float32`, but not
  ## an `unchained` unit (does not support `pow`)
  FloatLikeSupportsPow* = concept x, type T
    x.toFloat is float
    pow(x, 1.0) is T

  ## A concept for types that are float like and whose application of
  ## `^` with a RT integer value compiles and yields the same type. E.g. `float32`,
  ## but not an `unchained` unit (needs a `static int` to know the resulting
  ## type at CT)
  FloatLikeSupportsIntegerPowRT* = concept x, type T
    x.toFloat is float
    (x ^ 1) is T

  IdType = uint64 # "`uint64` should be enough for everyone... uhhh"

  DerivKey[T] = tuple[val, uncer: T, tag: IdType]

  Derivatives[T] = OrderedTable[DerivKey[T], T]

  ## *NOTE*: When dealing with `unchained` units, the derivative returned has the
  ## wrong units! Instead of `T` it should be the type of `∂a(x)/∂x`, but
  ## the code gets too complicated for the time being, as it would require use to
  ## store different types for the actual measurement and the gradient, turning
  ## `Measurement` into a double generic.
when (NimMajor, NimMinor, NimPatch) < (1, 7, 0):
  type
    Measurement*[T] = object
      val: T
      uncer: T
      id: IdType
      der: Derivatives[T] # map of the derivatives
else:
  type
    Measurement*[T: FloatLike] = object
      val*: T
      uncer*: T
      id: IdType
      der: Derivatives[T] # map of the derivatives


func value*[T: FloatLike](m: Measurement[T]): T {.inline.} = m.val
func uncertainty*[T: FloatLike](m: Measurement[T]): T {.inline.} = m.uncer
func error*[T: FloatLike](m: Measurement[T]): T {.inline.} = m.uncer

import hashes
proc hash*[T: FloatLike](key: DerivKey[T]): Hash =
  result = result !& hash(key.val)
  result = result !& hash(key.uncer)
  result = result !& hash(key.tag)
  result = !$result

when (NimMajor, NimMinor, NimPatch) < (1, 6, 0):
  ## add fallback `almostEqual`, backport from Nim std
  import fenv
  func almostEqual*[T: SomeFloat](x, y: T; unitsInLastPlace: Natural = 4): bool {.inline.} =
    if x == y:
      # short circuit exact equality -- needed to catch two infinities of
      # the same sign. And perhaps speeds things up a bit sometimes.
      return true
    let diff = abs(x - y)
    result = diff <= epsilon(T) * abs(x + y) * T(unitsInLastPlace) or
        diff < minimumPositiveValue(T)

  proc c_copysign(x, y: cdouble): cdouble {.importc: "copysign", header: "<math.h>".}
  func copySign*[T: SomeFloat](x, y: T): T {.inline.} =
    ## Returns a value with the magnitude of `x` and the sign of `y`;
    ## this works even if x or y are NaN, infinity or zero, all of which can carry a sign.
    result = c_copysign(x, y)

proc toFloat*[T: SomeFloat](x: T): float = x.float # float is biggest type, so this should be fine
proc toFloat*[T: SomeInteger](x: T): float = x.float # float is biggest type, so this should be fine

proc `==`*[T: FloatLike](k1, k2: DerivKey[T]): bool =
  if k1.tag == k2.tag and almostEqual(k1.val.float, k2.val.float) and
     almostEqual(k1.uncer.float, k2.uncer.float):
    result = true
  else:
    result = false

## Comparison operators. For two Measurements, only equal if that holds for
## value and uncertainty.
## TODO: should we use `almostEqual` here? Or is that too imprecise for the purpose?
proc `==`*[T: FloatLike](m1, m2: Measurement[T]): bool =
  ## comparison of two measurements does not need to take into account the
  ## type, as we require same type anyway. Hence use `toFloat` to compare as float
  mixin toFloat
  result = almostEqual(m1.val.toFloat, m2.val.toFloat) and
           almostEqual(m1.uncer.toFloat, m2.uncer.toFloat)

proc `==`*[T: FloatLike](m: Measurement[T], x: T): bool =
  result = almostEqual(m.val, x) and m.uncer == 0.0

proc `==`*[T: FloatLike](x: T, m: Measurement[T]): bool =
  result = almostEqual(m.val, x) and m.uncer == 0.0

proc isInf*[T: FloatLike](m: Measurement[T]): bool = m.val == Inf
proc isNegInf*[T: FloatLike](m: Measurement[T]): bool = m.val == -Inf

proc initDerivatives[T](): Derivatives[T] = initOrderedTable[DerivKey[T], T]()

proc changeType[T; U](tab: Derivatives[T], dtype: typedesc[U]): Derivatives[U] =
  result = initDerivatives[U]()
  for key, val in pairs(tab):
    let nkey = (val: key.val.U, uncer: key.uncer.U, tag: key.tag)
    result[nkey] = val.U

var idCounter = 0'u64
proc initMeasurement[T](val, uncer: T, der = initDerivatives[T](),
                        id = uint64.high): Measurement[T] =
  var id = if id == uint64.high: idCounter else: id # use id only if it's manually entered
  inc idCounter
  # create new `Derivatives` containing element itself if no elements in `der`
  let der = if der.len == 0: { (val: val, uncer: uncer, tag: id) : T(1) }.toOrderedTable
            else: der
  result = Measurement[T](val: val, uncer: uncer, id: id, der: der)

proc procRes[T](res: T, grad: T, arg: Measurement[T]): Measurement[T] =
  var der = initDerivatives[T]()
  for key, val in pairs(arg.der):
    if key.uncer != T(0.0):
      der[key] = (grad * val).T # force to `T` to avoid `unchained` errors
  # σ is 0 if input's σ is zero, else it's ∂(f(x))/∂x * Δx =^= grad * arg.uncer
  let σ = if arg.uncer == T(0.0): T(0.0) else: abs((grad * arg.uncer).T) # force to `T` to avoid `unchained` errors
  result = initMeasurement[T](res, σ, der, 0'u64)

proc derivative[T](a: Measurement[T], key: DerivKey[T]): T =
  ## See note about derivative in definition of `Measurement` type
  result = if key in a.der: a.der[key] else: T(0.0)

## A "type safe" solution required for `unchained` could be achieved with something like this.
## However, it shows the problem of getting a squared type as the `resder` argument
when false:
  proc procRes[T](res: T, grad: openArray[T], args: openArray[Measurement[T]]): Measurement[T] =
    assert grad.len == args.len, "`grad` and `args` must have the same number of elements!"
    var resder = initDerivatives[T]()
    #type SqrTyp = typeof(T(1) * T(1))
    var err = T(0) * T(0)
    for y in args:
      for key, val in y.der:
        if key notin resder:
          let σ_x = key.uncer
          if σ_x != T(0.0):
            var ∂G_∂x = T(0) * T(0)#: SqrTyp
            # iterate over all arguments of the function
            for (i, x) in pairs(args):
              # calc derivative of G w.r.t. current indep. var.
              let ∂a_∂x = derivative(x, key)
              if ∂a_∂x != T(0.0): # skip values with 0 partial deriv.
                ∂G_∂x = ∂G_∂x + (grad[i] * ∂a_∂x) # convert product back to avoid `unchained` error

            if ∂G_∂x != T(0) * T(0): #SqrTyp(0.0):
              # add (σ_x•∂G/∂x)² to total uncertainty (squared), but only if deriv != 0
              resder[key] = ∂G_∂x
              err = err + ((σ_x * ∂G_∂x) * (σ_x * ∂G_∂x)) # convert product back to avoid `unchained` error
    result = initMeasurement[T](res, sqrt(err), resder, 0'u64)

proc procRes[T](res: T, grad: openArray[T], args: openArray[Measurement[T]]): Measurement[T] =
  ## Note: In the body of this we perform rather funny back and forth conversions between `T`
  ## and sometimes float. This is because in `unchained` we generate new units for math operations.
  ## While this is fine in principle, it's too difficult to satisfy `unchained` CT safety logic
  ## here (or similar units `T` that change units on arithmetic). For example in theory the
  ## `Derivatives` would not have to store `T`, but rathe `T²`.
  ## This makes the implementation way to complex for no gain.
  assert grad.len == args.len, "`grad` and `args` must have the same number of elements!"
  var resder = initDerivatives[T]()
  var err: T
  for y in args:
    for key, val in pairs(y.der):
      if key notin resder:
        let σ_x = key.uncer
        if σ_x != T(0.0):
          var ∂G_∂x: T
          # iterate over all arguments of the function
          for (i, x) in pairs(args):
            # calc derivative of G w.r.t. current indep. var.
            let ∂a_∂x = derivative(x, key)
            if ∂a_∂x != T(0): # skip values with 0 partial deriv.
              ∂G_∂x = (∂G_∂x + (grad[i] * ∂a_∂x).T).T # convert product back to avoid `unchained` error
                                                  # technically this is a product after all
          if ∂G_∂x != T(0):
            # add (σ_x•∂G/∂x)² to total uncertainty (squared), but only if deriv != 0
            resder[key] = ∂G_∂x
            err = (err + ((σ_x * ∂G_∂x) * (σ_x * ∂G_∂x)).T).T # convert product back to avoid `unchained` error
  result = initMeasurement[T](res,
                              T(sqrt(err.float)), # convert to float and back to T to satisfy `unchained`
                              resder, 0'u64)

proc `±`*[T: FloatLike](val, uncer: T): Measurement[T] =
  result = initMeasurement[T](val, uncer)

proc `+-`*[T: FloatLike](val, uncer: T): Measurement[T] = val ± uncer
  ## `noUnicode`-mode makes output like this which users may want to re-use as
  ## input without edit.  Nim 2/`--experimental:unicodeOperators` allows it.

## Do we want the following? It forces any `FloatLike` to generate a `Measurement[float]`!
proc `±`*[T: FloatLike](val, uncer: T{lit}): Measurement[float] =
  result = initMeasurement[float](val.float, uncer.float)

proc measurement*[T: FloatLike](value, uncertainty: T): Measurement[T] =
  result = initMeasurement[T](value, uncertainty)

proc measurement*[T: FloatLike](val, uncer: T{lit}): Measurement[float] =
  result = initMeasurement[float](val.float, uncer.float)

when defined(useCligen):
  import cligen/strUt
  proc pretty*[T: FloatLike](m: Measurement[T], precision: int, merge = false): string =
    ## Uses `cligen` to format uncertainties. If `merge = false` uncertainties are
    ## printed in the form `A ± B` or `(A ± B)e-C`. If `merge = true` the uncertainty `B`
    ## is given as `A(B)e-C`.
    ##
    ## `merge = true` is useful as it allows to directly paste the output into `siunitx` in
    ## LaTeX.
    if merge:
      result = fmtUncertainMerged(m.val.float, m.uncer.float, sigDigs = precision)
    else:
      result = fmtUncertain(m.val.float, m.uncer.float, sigDigs = precision)
    when not (T is float): result.add " " & $T
else:
  proc pretty*[T: FloatLike](m: Measurement[T], precision: int): string =
    let mval = m.val.float.formatBiggestFloat(precision = precision)
    let merr = m.uncer.float.formatBiggestFloat(precision = precision)
    when not defined(noUnicode):
      result = &"{mval} ± {merr}"
    else:
      result = &"{mval} +- {merr}"
    when not (T is float):
      result.add " " & $T

proc `$`*[T: FloatLike](m: Measurement[T]): string = pretty(m, 3)

template print*(arg: untyped): untyped =
  echo astToStr(arg), ": ", $arg

template genOverloadsPlusMinus(fn: untyped): untyped =
  proc `fn`*[T: FloatLike](x: T, m: Measurement[T]): Measurement[T] = result = procRes(fn(x, m.val), T(1), m)
  proc `fn`*[T: FloatLike](m: Measurement[T], x: T): Measurement[T] = result = procRes(fn(m.val, x), T(1), m)
  ## Overloads for literals that force same type as Measurement has
  proc `fn`*[T: FloatLike; U: FloatLike](x: T{lit}, m: Measurement[U]): Measurement[U] =
    result = procRes(fn(U(x), m.val), U(1), m)
  proc `fn`*[U: FloatLike; T: FloatLike](m: Measurement[U], x: T{lit}): Measurement[U] =
    result = procRes(fn(m.val, U(x)), U(1), m)

## propagation based plainly on operator overload?
proc `+`*[T: FloatLike](a, b: Measurement[T]): Measurement[T] =
  result = procRes(T(a.val + b.val), [T(1), T(1)], [a, b])
# generate helpers
genOverloadsPlusMinus(`+`)

## mutable assign/addition
proc `+=`*[T: FloatLike](a: var Measurement[T], b: Measurement[T]) =
  let tmp = a + b
  a = tmp

proc `-`*[T: FloatLike](a, b: Measurement[T]): Measurement[T] =
  result = procRes(a.val - b.val, [T(1), T(-1)], [a, b])
# generate helpers
genOverloadsPlusMinus(`-`)

## mutable assign/subtraction
proc `-=`*[T: FloatLike](a: var Measurement[T], b: Measurement[T]) =
  let tmp = a - b
  a = tmp


## Type conversion. TODO: make this more type safe, funnily enough
proc to*[T: FloatLike; U](m: Measurement[T], dtype: typedesc[U]): Measurement[U] =
  when T is U:
    result = m
  else:
    result = initMeasurement[U](m.val.U, m.uncer.U, m.der.changeType(dtype), m.id)

# unary minus
proc `-`*[T: FloatLike](m: Measurement[T]): Measurement[T] = procRes(-m.val, T(-1), m)

proc `*`*[T: FloatLike; U: FloatLike](a: Measurement[T], b: Measurement[U]): auto =
  let val = a.val * b.val # TODO: the same logic should apply for the other cases.
                          # the equivalent line ensures only valid types are mixed.
  type resType = typeof(val)
  result = procRes(val, [b.val.resType, a.val.resType], [a.to(resType), b.to(resType)])

## The following two dirtry (!) templates help us with the assignment of the result
## procedure calls, taking care of deducing the types and performing the conversions
## so that we don't need to manually convert them all the time.
## `assign1` is used for the call to `procRes` with a single argument and
## `assign2` for the `procRes` with two arguments.
template assign1(valArg, arg, m: untyped): untyped {.dirty.} =
  let val = valArg
  type resType = typeof(val)
  result = procRes(val, arg.resType, m.to(resType))

template assign2(valArg, arg1, arg2, m1, m2: untyped): untyped {.dirty.} =
  let val = valArg
  type resType = typeof(val)
  result = procRes(val, [arg1.resType, arg2.resType], [m1.to(resType), m2.to(resType)])

# helper overloads
proc `*`*[T: FloatLike; U: FloatLike](x: T, m: Measurement[U]): auto =
  assign1(x * m.val, x, m)

proc `*`*[T: FloatLike; U: FloatLike](m: Measurement[U], x: T): auto =
  assign1(x * m.val, x, m)

## Overloads for literals that force same type as Measurement has
proc `*`*[T: FloatLike; U: FloatLike](x: T{lit}, m: Measurement[U]): Measurement[U] =
  result = procRes(U(U(x) * m.val), U(x), m)
proc `*`*[U: FloatLike; T: FloatLike](m: Measurement[U], x: T{lit}): Measurement[U] =
  result = procRes(U(m.val * U(x)), U(x), m)

## mutable assign/multiplication
proc `*=`*[T: FloatLike](a: var Measurement[T], b: Measurement[T]) =
  let tmp = a * b
  a = tmp

proc `*=`*[T: FloatLike](a: var Measurement[T], b: T) =
  let tmp = a * b
  a = tmp


proc `/`*[T: FloatLike; U: FloatLike](a: Measurement[T], b: Measurement[U]): auto =
  let oneOverB = 1.0 / b.val
  assign2(a.val / b.val, oneOverB, -a.val * (oneOverB * oneOverB), a, b)

# helper overloads
proc `/`*[T: FloatLike; U: FloatLike](x: T, m: Measurement[U]): auto =
  assign1(x / m.val, -x / (m.val * m.val), m)

proc `/`*[T: FloatLike; U: FloatLike](m: Measurement[T], x: U): auto =
  assign1(m.val / x, 1.0 / x, m)

## Overloads for literals that force same type as Measurement has
proc `/`*[T: FloatLike; U: FloatLike](x: T{lit}, m: Measurement[U]): auto =
  when T is SomeInteger:
    let x = U(x)
  type A = typeof( x / m.val )
  let arg: A = A(x / m.val)
  let grad: A =  A(-x / (m.val * m.val))
  result = procRes(arg, grad, m.to(A))
proc `/`*[U: FloatLike; T: FloatLike](m: Measurement[U], x: T{lit}): Measurement[U] =
  when T is SomeInteger:
    let x = U(x)
  result = procRes(U(m.val / x), 1.0 / x, m)

# Power `^`
## NOTE: Using any of the following exponentiation functions is dangerous. The Nim parser
## might include an additional prefix into the argument to `^` instead of keeping it as a,
## well, prefix. So `-(x - μ)^2` is parsed as `(-(x - μ))^2`!

from measuremancer / math_utils import power
## -> for static exponents
proc `**`*[T: FloatLike](m: Measurement[T], p: static SomeInteger): auto =
  # Get the resulting type
  type U = typeof(power(m.val, p))
  # and convert what is necessary.
  result = procRes(U(power(m.val, p)), U(p.float * power(m.val, (p - 1))), m.to(U))
proc `^`*[T: FloatLike](m: Measurement[T], p: static SomeInteger): auto =
  ## NOTE: If you import `unchained`, this version is not actually used, as `unchained`
  ## defines a macro `^` for static integer exponents, which simply rewrites
  ## the AST to an infix (or trivial) node!
  m ** p

## -> for RT natural exponents
proc `**`*[T: FloatLikeSupportsIntegerPowRT](m: Measurement[T], p: SomeInteger): Measurement[T] =
  result = procRes(power(m.val, p), p.float * power(m.val, (p - 1)), m)
proc `^`*[T: FloatLikeSupportsIntegerPowRT](m: Measurement[T], p: SomeInteger): Measurement[T] =
  m ** p

## -> for explicit float exponents
proc `**`*[T: FloatLikeSupportsPow; U: FloatLikeNotInt](
  m: Measurement[T], p: U): Measurement[T] =
    result = procRes(pow(m.val, p), p * pow(m.val, (p - 1.0)), m)
proc `^`*[T: FloatLikeSupportsPow; U: FloatLikeNotInt](
  m: Measurement[T], p: U): Measurement[T] =
    m ** p

proc `**`*[T: FloatLike](a, b: Measurement[T]): Measurement[T] =
  let powV = pow(a.val, b.val)
  result = procRes(powV,
                   [pow(a.val, b.val - 1.0) * b.val,
                    powV * ln(a.val)],
                   [a, b])
proc `^`*[T: FloatLike](a, b: Measurement[T]): Measurement[T] = a ** b

# TODO: need to add a lot more primitive functions

## Now generate single argument (mostly trigonometric) functions
## They all follow the following structure
## ```
##   proc exp*(m: Measurement): Measurement =
##     result = procRes(exp(m.val), exp(m.val), m)
## ```
## where the derivatives are taken from the 'lookup table' below.

import std / macros
template genCall(fn, letStmt, x, m, deriv: untyped): untyped =
  proc `fn`*[T](m: Measurement[T]): auto =
    ## letStmt contains the definition of x, `let x = m.val`
    letStmt
    ## `deriv` is handed as an expression containing `x`
    type U = typeof(fn(x))
    result = procRes(fn(x), deriv, m.to(U))

macro defineSupportedFunctions(body: untyped): untyped =
  result = newStmtList()
  for fn in body:
    doAssert fn.kind == nnkInfix and fn[0].strVal == "->"
    let fnName = fn[1].strVal
    let fnId = ident(fnName)
    # generate the code
    let deriv = fn[2]
    let xId = ident"x"
    let mId = ident"m"
    let xLet = quote do:
      let `xId` = `mId`.val
    let ast = getAst(genCall(fnId, xLet, xId, mId, deriv))
    result.add ast

## NOTE: some of the following functions are not implemented in Nim atm, hence
## they are commented out.
## This is based on the same in `astgrad`, which itself took the map from
## somewhere else. Need to look up where that was, oops.
defineSupportedFunctions:
  sqrt        ->  1.0 / 2.0 / sqrt(x)
  cbrt        ->  1.0 / 3.0 / (cbrt(x)^2.0)
  abs2        ->  1.0 * 2.0 * x
  inv         -> -1.0 * abs2(inv(x))
  log         ->  1.0 / x
  log10       ->  1.0 / x / log(10)
  log2        ->  1.0 / x / log(2.0)
  log1p       ->  1.0 / (x + 1.0)
  exp         ->  exp(x)
  exp2        ->  log(2.0) * exp2(x)
  expm1       ->  exp(x)
  sin         ->  cos(x)
  cos         -> -sin(x)
  tan         ->  (1.0 + (tan(x)^2))
  sec         ->  sec(x) * tan(x)
  csc         -> -csc(x) * cot(x)
  cot         -> -(1.0 + (cot(x)^2))
  #sind        ->  Pi / 180.0 * cosd(x)
  #cosd        -> -Pi / 180.0 * sind(x)
  #tand        ->  Pi / 180.0 * (1.0 + (tand(x)^2))
  #secd        ->  Pi / 180.0 * secd(x) * tand(x)
  #cscd        -> -Pi / 180.0 * cscd(x) * cotd(x)
  #cotd        -> -Pi / 180.0 * (1.0 + (cotd(x)^2))
  arcsin      ->  1.0 / sqrt(1.0 - (x^2))
  arccos      -> -1.0 / sqrt(1.0 - (x^2))
  arctan      ->  1.0 / (1.0 + (x^2))
  arcsec      ->  1.0 / abs(x) / sqrt(x^2 - 1.0)
  arccsc      -> -1.0 / abs(x) / sqrt(x^2 - 1.0)
  arccot      -> -1.0 / (1.0 + (x^2))
  #arcsind     ->  180.0 / Pi / sqrt(1.0 - (x^2))
  #arccosd     -> -180.0 / Pi / sqrt(1.0 - (x^2))
  #arctand     ->  180.0 / Pi / (1.0 + (x^2))
  #arcsecd     ->  180.0 / Pi / abs(x) / sqrt(x^2 - 1.0)
  #arccscd     -> -180.0 / Pi / abs(x) / sqrt(x^2 - 1.0)
  #arccotd     -> -180.0 / Pi / (1.0 + (x^2))
  sinh        ->  cosh(x)
  cosh        ->  sinh(x)
  tanh        ->  sech(x)^2
  sech        -> -tanh(x) * sech(x)
  csch        -> -coth(x) * csch(x)
  coth        -> -(csch(x)^2)
  arcsinh     ->  1.0 / sqrt(x^2 + 1.0)
  arccosh     ->  1.0 / sqrt(x^2 - 1.0)
  arctanh     ->  1.0 / (1.0 - (x^2))
  arcsech     -> -1.0 / x / sqrt(1.0 - (x^2))
  arccsch     -> -1.0 / abs(x) / sqrt(1.0 + (x^2))
  arccoth     ->  1.0 / (1.0 - (x^2))
  deg2rad     ->  Pi / 180.0
  rad2deg     ->  180.0 / Pi
  erf         ->  2.0 * exp(-x*x) / sqrt(Pi)
  erfinv      ->  0.5 * sqrt(Pi) * exp(erfinv(x) * erfinv(x))
  erfc        -> -2.0 * exp(-x*x) / sqrt(Pi)
  erfcinv     -> -0.5 * sqrt(Pi) * exp(erfcinv(x) * erfcinv(x))
  erfi        ->  2.0 * exp(x*x) / sqrt(Pi)
  #gamma       ->  digamma(x) * gamma(x)
  #lgamma      ->  digamma(x)
  #digamma     ->  trigamma(x)
  #invdigamma  ->  inv(trigamma(invdigamma(x)))
  #trigamma    ->  polygamma(2.0 x)
  #airyai      ->  airyaiprime(x)
  #airybi      ->  airybiprime(x)
  #airyaiprime ->  x * airyai(x)
  #airybiprime ->  x * airybi(x)
  #besselj0    -> -besselj1(x)
  #besselj1    ->  (besselj0(x) - besselj(2.0, x)) / 2.0
  #bessely0    -> -bessely1(x)
  #bessely1    ->  (bessely0(x) - bessely(2.0, x)) / 2.0
  #erfcx       ->  (2.0 * x * erfcx(x) - 2.0 / sqrt(Pi))
  #dawson      ->  (1.0 - 2.0 * x * dawson(x))

func signbit*[T: FloatLike](m: Measurement[T]): bool = m.val.signbit

proc copysign*[T: FloatLike](a, b: Measurement[T]): Measurement[T] =
  result = if signbit(a) != signbit(b): -a else: a

proc copysign*[T: FloatLike; U: FloatLike](a: Measurement[T], b: U): Measurement[T] =
  result = if signbit(a) != signbit(b): -a else: a
proc copysign*[T: FloatLike; U: FloatLike](a: U, b: Measurement[T]): U =
  result = if signbit(a) != signbit(b): -a else: a

proc abs*[T: FloatLike](m: Measurement[T]): Measurement[T] =
  let mval = m.val
  result = procRes(abs(mval), copysign(mval, 1.T), m)

template comp(fn: untyped): untyped =
  proc `fn`*[T: FloatLike](a, b: Measurement[T]):    bool = result = fn(a.val, b.val)
  proc `fn`*[T: FloatLike](a: Measurement[T], b: T): bool = result = fn(a.val, b)
  proc `fn`*[T: FloatLike](a: T, b: Measurement[T]): bool = result = fn(a, b.val)
  proc `fn`*[T: FloatLike; U: FloatLike](a: Measurement[T], b: U{lit}): bool = result = fn(a.val, T(b))
  proc `fn`*[T: FloatLike; U: FloatLike](a: U{lit}, b: Measurement[T]): bool = result = fn(T(a), b.val)

comp(`<`)
comp(`<=`)

when isMainModule:
  proc foo[T](x: T, μ, σ: float): T =
    let val = 2 * σ * σ
    result = exp(- (x*x) / val ) #/ (2 * σ^2))

  let x1 = 5.0 ± 1.0
  let x2 = 2.5 ± 1.5

  let y = 2 ± 2
  print y
  print x1

  print x1 + x2
  print x1 - x2
  print x1 * x2
  print x1 / x2
  print exp(x1)
  print exp(2.0 * x1)
  print x1 * x1
  print x1 * x1 * x1
  print exp(x1 * x1)

  print 2.0 * x1
  print 2 / x1
  print x1 / 2


  when false:
    import unchained
    let k1 = 5.0.keV ± 1.0.keV
    let k2 = 2.5.keV ± 1.5.keV
    print k1 + k2

    print foo(x1, 33.333, 100.0)
    print exp(- x1*x1 / (2 * 100.0 * 100.0))

  when false:
    import arraymancer
    let t = [x1, x2].toTensor.map_inline(x * x)

    for x in t:
      echo x


  #
  print x1^2
  #
  #let x3 = 1.0 ± 1.0
  #let x4 = 1.0 ± 1.0
  #echo x3 / x3
  #echo x3 - x4

  #let nofloatbutconv = 1 ± 1
  #echo nofloatbutconv

  #let nofloatnoconv = "a" ± "b"
  #echo nofloatnoconv
