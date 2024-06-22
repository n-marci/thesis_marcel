set TEXINPUTS=.\themefau

latex Masterarbeit.tex
bitex Masterarbeit
latex Masterarbeit.tex
latex Masterarbeit.tex

del *.aux *.bbl *.bcf *.blg *.nav *.out *.snm *.log *.toc *.vrb
exit
