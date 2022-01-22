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

  IdType = uint64 # "`uint64` should be enough for everyone... uhhh"

  DerivKey[T] = tuple[val, uncer: T, tag: IdType]

  Derivatives[T] = OrderedTable[DerivKey[T], T]

  Measurement[T: FloatLike] = object #{.requiresInit.} = object
    val: T
    uncer: T
    id: IdType
    der: Derivatives[T] # map of the derivatives

import hashes
proc hash*[T: FloatLike](key: DerivKey[T]): Hash =
  result = result !& hash(key.val)
  result = result !& hash(key.uncer)
  result = result !& hash(key.tag)
  result = !$result

proc `==`*[T: FloatLike](k1, k2: DerivKey[T]): bool =
  if k1.tag == k2.tag and almostEqual(k1.val.float, k2.val.float) and
     almostEqual(k1.uncer.float, k2.uncer.float):
    result = true
  else:
    result = false

proc toFloat*[T: SomeFloat](x: T): float = x.float # float is biggest type, so this should be fine

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
  result = initMeasurement[T](res, σ, der, 0'u32)

proc derivative[T](a: Measurement[T], key: DerivKey[T]): T =
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
    result = initMeasurement[T](res, sqrt(err), resder, 0'u32)

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
              ∂G_∂x = ∂G_∂x + (grad[i] * ∂a_∂x).T # convert product back to avoid `unchained` error
                                                  # technically this is a product after all
          if ∂G_∂x != T(0):
            # add (σ_x•∂G/∂x)² to total uncertainty (squared), but only if deriv != 0
            resder[key] = ∂G_∂x
            err = err + ((σ_x * ∂G_∂x) * (σ_x * ∂G_∂x)).T # convert product back to avoid `unchained` error
  result = initMeasurement[T](res,
                              T(sqrt(err.float)), # convert to float and back to T to satisfy `unchained`
                              resder, 0'u32)

proc `±`*[T: FloatLike](val, uncer: T): Measurement[T] =
  result = initMeasurement[T](val, uncer)

## Do we want the following? It forces any `FloatLike` to generate a `Measurement[float]`!
proc `±`*[T: FloatLike](val, uncer: T{lit}): Measurement[float] =
  result = initMeasurement[float](val.float, uncer.float)

proc measurement*[T: FloatLike](value, uncertainty: T): Measurement[T] =
  result = initMeasurement[float](value, uncertainty)

proc pretty*[T: FloatLike](m: Measurement[T], precision: int): string =
  let mval = m.val.float.formatBiggestFloat(precision = precision)
  let merr = m.uncer.float.formatBiggestFloat(precision = precision)
  result = &"{mval} ± {merr}"
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
    result = procRes(fn(x, m.val), U(1), m)
  proc `fn`*[U: FloatLike; T: FloatLike](m: Measurement[U], x: T{lit}): Measurement[U] =
    result = procRes(fn(m.val, x), U(1), m)

## propagation based plainly on operator overload?
proc `+`*[T: FloatLike](a, b: Measurement[T]): Measurement[T] =
  result = procRes(T(a.val + b.val), [T(1), T(1)], [a, b])
# generate helpers
genOverloadsPlusMinus(`+`)

proc `-`*[T: FloatLike](a, b: Measurement[T]): Measurement[T] =
  result = procRes(a.val - b.val, [T(1), T(-1)], [a, b])
# generate helpers
genOverloadsPlusMinus(`-`)

## Type conversion. TODO: make this more type safe, funnily enough
proc to[T: FloatLike; U](m: Measurement[T], dtype: typedesc[U]): Measurement[U] =
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

proc `*`*[T: FloatLike](m: Measurement[T], x: T): Measurement[T] =
  assign1(m.val * x, x, m)

## Overloads for literals that force same type as Measurement has
proc `*`*[T: FloatLike; U: FloatLike](x: T{lit}, m: Measurement[U]): Measurement[U] =
  result = procRes(U(x) * m.val, U(x), m)
proc `*`*[U: FloatLike; T: FloatLike](m: Measurement[U], x: T{lit}): Measurement[U] =
  result = procRes(m.val * U(x), U(x), m)


proc `/`*[T: FloatLike; U: FloatLike](a: Measurement[T], b: Measurement[U]): auto =
  let oneOverB = 1.0 / b.val
  assign2(a.val / b.val, oneOverB, -a.val * (oneOverB * oneOverB), a, b)

# helper overloads
proc `/`*[T: FloatLike; U: FloatLike](x: T, m: Measurement[U]): auto =
  assign1(x / m.val, -x / (m.val * m.val), m)

proc `/`*[T: FloatLike; U: FloatLike](m: Measurement[T], x: U): auto =
  assign1(m.val / x, 1.0 / x, m)

## Overloads for literals that force same type as Measurement has
proc `/`*[T: FloatLike; U: FloatLike](x: T{lit}, m: Measurement[U]): Measurement[U] =
  result = procRes(U(x) / m.val, U(-x) / (m.val * m.val), m)
proc `/`*[U: FloatLike; T: FloatLike](m: Measurement[U], x: T{lit}): Measurement[U] =
  result = procRes(m.val / U(x), 1.0 / U(x), m)

# Power `^`
## NOTE: Using any of the following exponentiation functions is dangerous. The Nim parser
## might include an additional prefix into the argument to `^` instead of keeping it as a,
## well, prefix. So `-(x - μ)^2` is parsed as `(-(x - μ))^2`!
proc `**`*[T: FloatLike](m: Measurement[T], p: Natural): Measurement[T] =
  result = procRes(m.val ^ p, p.float * m.val ^ (p - 1), m)
proc `^`*[T: FloatLike](m: Measurement[T], p: Natural): Measurement[T] = m ** p

proc `**`*[T: FloatLike](m: Measurement[T], p: FloatLike): Measurement[T] =
  result = procRes(pow(m.val, p), p * pow(m.val, (p - 1.0)), m)
proc `^`*[T: FloatLike](m: Measurement[T], p: FloatLike): Measurement[T] = m ** p


proc `**`*[T: FloatLike](a, b: Measurement[T]): Measurement[T] =
  let powV = pow(a.val, b.val)
  result = procRes(powV,
                   [pow(a.val, b.val - 1.0),
                    powV * log(a.val)],
                   [a, b])
proc `^`*[T: FloatLike](a, b: Measurement[T]): Measurement[T] = a ** b

# TODO: need to add a lot more primitive functions
# Primitives
proc exp*(m: Measurement): Measurement =
  result = procRes(exp(m.val), exp(m.val), m)



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
