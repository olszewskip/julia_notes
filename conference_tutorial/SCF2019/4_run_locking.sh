#!/bin/bash

export JULIA_NUM_THREADS=1
julia 4_locking.jl

export JULIA_NUM_THREADS=$(nproc)
julia 4_locking.jl

