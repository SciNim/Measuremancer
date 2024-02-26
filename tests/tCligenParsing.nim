import std / random
## This is just a short program to generate numbers (value ± uncertainty),
## format them with `cligen` and verify that the formatting is what we expect, given
## the uncertainties and the required precision.
import ../measuremancer
import std / strscans

template checkError(typ, body, msg: untyped): untyped =
  if not body:
    inc typ
    echo msg

var
  errUncert = 0
  errDigits = 0
  errExp = 0
  errDigit = 0

# How many digits the error may be smaller than the value in OOM
const N_smaller = 6
const Total = 1_000_000
for i in 0 ..< Total:
  if i mod 50_000 == 0:
    echo "i = ", i, " / ", Total
  # 1. sample exponents
  let expVal = rand(-300 .. 300)
  # 1b. exponent up to N digits smaller
  let expErrDiff = rand(0 .. N_smaller)
  # 2. sample a value and an eror
  let valNoExp = rand(0.1 .. 1.0)
  let errNoExp = rand(0.1 .. 1.0)
  # 3. combine value and error
  let val = valNoExp * 10 * pow(10.0, (expVal).float)
  let err = errNoExp * 10 * pow(10.0, (expVal - expErrDiff).float)
  # 4. -> make a 'Measurement'
  let meas = val ± err
  # 5. sample a precision
  let prec = rand(1 .. 16)
  let resStr = pretty(meas, precision = prec, merge = true)
  # 6. verify the produced string
  #echo val, " ± ", err, " = ", resStr, " prec = ", prec
  # Note that we parse digits after `.` and precision as strings. Might exceed `int.high` and we don't
  # want to add an additional int ⇔ string conversion
  let (success, digit, numbers, uncert, exp) = scanTuple(resStr, "$i.$+($+)e$i")
  if success:
    # 1. verify precision in parenthesis. Must be number of digits given by `prec`
    checkError errUncert, (uncert).len == prec, ("Need uncertainty: " & $prec & " got: " & $uncert)
    # 2. determine difference in exponent after possible rounding
    var expErrReal = expErrDiff
    # 2a. check if rounding of error affects difference
    let errNoExpRound = errNoExp.round(prec) # Rounding of error can affect rounding of value
    if errNoExpRound == 1.0: # increased by OOM, adjust difference in exponent
      expErrReal -= 1
    # 2b. round according to (possibly) updated exp difference...
    let valNoExpRound = valNoExp.round(expErrReal + prec) # use final exponent difference and precision to round
    let valueRound = expErrReal + prec # store this value to round integer by
    # 2c. and possibly adjust exponent difference
    if valNoExpRound == 1.0: # increased by OOM, adjust exp difference
      expErrReal += 1
    #echo valNoExpRound, " vs ", valNoExp, " and ", errNoExpRound, " vs ", errNoExp
    # 2d. `-1` for one significant digit in front of `.`
    expErrReal -= 1
    let numDigits = min(expErrReal + prec, 18) # maximum precision: 18
    #echo "expErrDiff ", expErrDiff, " expErrReal ", expErrReal, " prec ", prec, " numDigits ", numDigits, " numbers: ", numbers.len
    #echo "REAL EXP: ", expErrDiff + prec - 1
    checkError errDigits, numDigits == numbers.len, ("Need digits: " & $numDigits & " got: " & $numbers.len)
    # 3. verify exponent
    let expValReal = if valNoExpRound == 1.0: expVal + 1 else: expVal # if value rounded up!
    checkError errExp, exp == expValReal, ("Need exponent: " & $expValReal & " got: " & $exp)
    # 4. verify digit in front of `.`
    var intDigit = ((valNoExp * 10.0).round(valueRound - 1)).trunc.int
    if intDigit == 10:
      intDigit = 1
    checkError errDigit, intDigit == digit, "Need integer: " & $intDigit & " got: " & $digit

  if not success: ## Cases without `e` notation
    # HANDLE THES
    #echo success, " res ", resStr
    discard

echo "Found errors:"
echo "\t Uncertainty = ", errUncert
echo "\t Digits      = ", errDigits
echo "\t Exponent    = ", errExp
echo "\t Digit       = ", errDigit

doAssert errUncert == 0 and errDigits == 0 and errExp == 0 and errDigit == 0
