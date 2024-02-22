# Package

version       = "0.2.6"
author        = "Vindaar"
description   = "A library to handle measurement uncertainties"
license       = "MIT"

# Dependencies

requires "nim >= 1.4.0"

task testDeps, "Installs dependencies for tests":
  exec "nimble install -y unchained"

task test, "Run standard tests":
  exec "nim c -r tests/tmeasuremancer.nim"
  exec "nim c -d:release -d:useCligen -r tests/tCligenParsing.nim"
