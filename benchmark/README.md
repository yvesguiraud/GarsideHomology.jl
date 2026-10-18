# Benchmarks

Measurement harnesses and protocols for Gauss.jl; they are not part of the loaded package.

Protocol: one fresh Julia session per measurement, at least three runs in random order,
minimum kept; machine, Julia version and commit written in every CSV line; results saved
incrementally. Timings are compared only on the same machine and Julia version.
