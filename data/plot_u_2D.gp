type = ARG1
dir = "data/"
data_file = system("ls -1v " . dir . "profile_2D.dat")

set terminal qt
unset key
set view map
set xrange [*:*] noreverse writeback
set x2range [*:*] noreverse writeback
set yrange [*:*] noreverse writeback
set y2range [*:*] noreverse writeback
set zrange [*:*] noreverse writeback
set cbrange [*:*] noreverse writeback
set rrange [*:*] noreverse writeback
set dgrid3d 128,128 splines
set palette @MATLAB
set xlabel 'x'
set ylabel 'y'
set zlabel 'u'
set border 3
set xrange [0.0:1.0]
set yrange [0.0:1.0]
set xtics 0.0, 0.1, 1.0
set ytics 0.0, 0.1, 1.0
set xtics nomirror
set ytics nomirror
set mxtics 5
set mytics 5
# set key spacing 1.3 font "Arial, 9"
# set key box lt -1 lw 1.0

splot data_file u 1:2:3 w pm3d t "SLDG" 

output_file = dir . 'profile_2D'
if ((type eq '') || (type eq 'png')) {
     term = 'pngcairo'
     ext = '.png'
     dpi = 600 
     width = 7.2 ## inch (variable)
     height = 5.4 ## inch (variable)
     pt2in = 0.01389 
     in2px = dpi
     ptscale = pt2in * in2px * 0.75
     round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
     wpx = round(width * in2px)
     hpx = round(height * in2px)
     output_file = output_file. ext
     set terminal term size wpx,hpx fontscale ptscale linewidth ptscale pointscale ptscale
}
if (type eq 'pdf') {
     term = 'pdfcairo'
     ext = '.pdf'
     output_file = output_file . ext
     set terminal term size 7.2in,5.4in fontscale 0.75 lw 1.5 pointscale 0.75
}
if (type eq 'eps') { 
     term = 'epscairo'
     ext = '.eps'
     output_file = output_file . ext
     set terminal term size 7.2in,5.4in fontscale 0.75 lw 1.5 pointscale 0.75
}
if (type eq 'tex') { 
     term = 'cairolatex'
     ext = '.tex' 
     output_file = output_file . ext
     set terminal term pdf size 7.2in,5.4in fontscale 0.75 lw 1.5 pointscale 0.75
}

set output output_file
replot
set output
set terminal qt

