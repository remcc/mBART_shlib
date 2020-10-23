The following should work in linux (tested in Ubuntu 18.04.5 LTS).
I have also got it to work easily on my MacBook Pro with the xcode installed.

Make sure you have the R package Rcpp installed. 

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
I am using section 2.5 of the Rcpp FAQ (frequently asked questions).
Dirk will scorn this, but is very easy.
