#!/usr/bin/env gnuplot
#
# This sample script plots a generic bunch of graphs in a single figure.
#


NCOLS=3
NROWS=6
set tmargin 0.5
set bmargin 3
set lmargin 7
set rmargin 0.5

set term post enh c eps size NCOLS,NROWS
set out 'multiplot.eps'
set multiplot layout NROWS,NCOLS


do for [icol=1:NCOLS] {
	do for [irow=1:NROWS] {
		unset key
		set title "icol=".icol." irow=".irow
        p[-1:1] sin(x*icol+irow)
    }
}
unset multiplot
