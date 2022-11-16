import std / macros

# Taken from `unchained`. This won't be exported from `measuremancer` though,
# rather we will `bind` in the scope in which we use this.
macro power*(x: typed, num: static int): untyped =
  ## general purpose power using `^` for integers, which works for any
  ## type by rewriting to a product of `*`.
  ##
  ## For the special cases of -1, 0, 1 we simply rewrite to the correct
  ## result. Negative powers are written as `1 / x^p`
  if num == 0:
    result = quote do:
      1.0
  elif num == 1:
    result = x
  elif num == -1:
    ## Assume that the type supports addition by a float!
    result = quote do:
      1.0 / `x`
  else:
    result = nnkInfix.newTree(ident"*")

    proc addInfix(n, x: NimNode, num: int) =
      var it = n
      if num > 0:
        it.add nnkInfix.newTree(ident"*")
        it[1].addInfix(x, num - 1)
      while it.len < 3:
        it.add x

    result.addInfix(x, abs(num) - 2)

    # invert if num is negative
    if num < -1:
      ## Assume that the type supports addition by a float!
      result = quote do:
        1.0 / (`result`)

import std/math
proc power*[T](x: T, num: SomeInteger): T =
  if num > 0: result = x ^ num
  elif num == 0: result = T(1)
  else: result = T(1) / (x ^ abs(num))
