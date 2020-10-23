The following may work in linux (tested in Ubuntu 18.04.5 LTS).

(i)
As indicated in RMakefile do this at the linux command line
$ export PKG_CXXFLAGS=`Rscript -e "Rcpp:::CxxFlags()"`

(ii) at the linux command line:
 $ make 

- you should see the files load-monbart.R and monbart.so

(iii)
 work through the example script test.R

Have a look at load-monbart.R.
The idea is that you can just copy this file to any directory you want to work in
and then use > source("load-monbart.R").

Note:
Instead of (i) and (ii) you can try:
$ make -f robMakefile
  but this assumes the paths on your machine are the same as Rob's current machine.
