# Package

version       = "0.1.1"
author        = "Vindaar"
description   = "A library to handle measurement uncertainties"
license       = "MIT"

# Dependencies

requires "nim >= 1.4.0"

task testCI, "Run standard tests in CI (installs unchained)":
  exec "nimble install -y unchained"
  exec "nim c -r tests/tmeasuremancer.nim"

task test, "Run standard tests":
  exec "nim c -r tests/tmeasuremancer.nim"
