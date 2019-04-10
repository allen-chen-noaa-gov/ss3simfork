all: plosaccepted

arxiv:
	pdflatex ss3sim-ms-arxiv
	pdflatex ss3sim-ms-arxiv
	pdflatex ss3sim-ms-arxiv

plosaccepted:
	pdflatex ss3sim-ms-plos.tex
	pdflatex ss3sim-ms-plos.tex
	pdflatex ss3sim-ms-plos.tex

plossubmit:
	pandoc ~/Dropbox/ms/ss3sim-response/ss3sim-response.md -o ss3sim-response.tex
	pdflatex ss3sim-ms
	#pdflatex ss3sim-ms
	rm ss3sim-response.tex
	cp ss3sim-ms.pdf ~/Dropbox/Public/

plos:
	pandoc -S --no-wrap --bibliography=ss3sim-ms.bib --natbib ss3sim-ms.md -o manuscript.tex
	perl -p -i -e "s/Fig. /Fig.~/g" manuscript.tex
	perl -p -i -e "s/citep/cite/g" manuscript.tex
	perl -p -i -e "s/citet/cite/g" manuscript.tex
	perl -p -i -e "s/e\.g\. /e\.g\.~/g" manuscript.tex
	perl -p -i -e "s/i\.e\. /i\.e\.~/g" manuscript.tex
	perl -p -i -e "s/vs. /vs.~/g" manuscript.tex
	perl -p -i -e "s/No. /No.~/g" manuscript.tex
	perl -p -i -e "s/\\\section/\\\section\*/g" manuscript.tex
	perl -p -i -e "s/\\\subsection/\\\subsection\*/g" manuscript.tex
	perl -p -i -e "s/Andre Punt/Andr\\\'{e} Punt/g" manuscript.tex
	pandoc ~/Dropbox/ms/ss3sim-response/ss3sim-response.md -o ss3sim-response.tex
	pdflatex ss3sim-ms
	bibtex ss3sim-ms
	pdflatex ss3sim-ms
	pdflatex ss3sim-ms
	rm ss3sim-response.tex
	cp ss3sim-ms.pdf ~/Dropbox/Public/

diff:
	#cp manuscript-plos.tex ms-temp.tex
	#perl -p -i -e "s/\\\textsf{R}/R/g" ms-temp.tex
	#perl -p -i -e "s/\\\textbf{ss3sim}/ss3sim/g" ms-temp.tex
	#perl -p -i -e "s/\\\textbf{r4ss}/r4ss/g" ms-temp.tex
	latexdiff ss3sim-first-submission.tex ss3sim-ms-for-diff.tex > ss3sim-ms-diff.tex
	pdflatex ss3sim-ms-diff
	#pdflatex ss3sim-ms-diff

letter:
	xelatex ss3sim-letter

bib:
	pdflatex ss3sim-ms
	bibtex ss3sim-ms
	pdflatex ss3sim-ms
	pdflatex ss3sim-ms

clean:
	rm *.log
	rm *.aux
	rm *.bbl
	rm *.blg

docx:
	pandoc -S --bibliography=ss3sim-ms.bib --csl=plos.csl --reference-docx=reference.docx ss3sim-ms.md -o ss3sim-ms.docx
	mv ss3sim-ms.docx ~/Dropbox/Public/

rft:
	latex2rtf ss3sim-ms.tex
	mv ss3sim-ms.rtf ~/Dropbox/Public/

md:
	pandoc -S --bibliography=ss3sim-ms.bib --csl=plos.csl ss3sim-ms.md -o ss3sim-ms-with-refs.md

cleanbib:
	bibtool ss3sim-ms.bib -s > ss3sim-ms2.bib
	mv ss3sim-ms2.bib ss3sim-ms.bib

vignette:
	./knit
	pandoc -S -N ss3sim-vignette.md -o ss3sim-vignette.pdf
	rm ss3sim-vignette.md
