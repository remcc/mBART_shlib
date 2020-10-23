## for this to work do this at the unix command line:
#   $ export PKG_CXXFLAGS=`Rscript -e "Rcpp:::CxxFlags()"`

all: monbart.so load-monbart.R

monbart.so: cmonbart.cpp funs.cpp tree.cpp bd.cpp
	R CMD SHLIB cmonbart.cpp funs.cpp tree.cpp bd.cpp -I/usr/local/lib/R/library/Rcpp/include
	mv cmonbart.so monbart.so

load-monbart.R:
	rm -f load-monbart.R
	echo "dyn.load(\"$(PWD)/monbart.so\")" > load-monbart.R
	echo "source(\"$(PWD)/rpmonbart.R\")" >> load-monbart.R
	echo "source(\"$(PWD)/monbart.R\")" >> load-monbart.R

clean: 
	rm -f monbart.so
	rm -f load-monbart.R
	rm -f *.o
