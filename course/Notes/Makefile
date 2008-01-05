# thelecture 0 = all, 1 = lecture1,...
all:
	-pdflatex  -jobname scicomp $(FLAGS) "\nonstopmode\newcount\thelecture\thelecture=0\relax\input{scicomp.tex}"
	-pdflatex  -jobname scicomp $(FLAGS) "\nonstopmode\newcount\thelecture\thelecture=0\relax\input{scicomp.tex}"



lecture1:
	-pdflatex  -jobname scicomp-lecture1 $(FLAGS) "\nonstopmode\newcount\thelecture\thelecture=1\relax\input{scicomp.tex}"
	-pdflatex  -jobname scicomp-lecture1 $(FLAGS) "\nonstopmode\newcount\thelecture\thelecture=1\relax\input{scicomp.tex}"
