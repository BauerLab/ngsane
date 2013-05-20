#%Module1.0#####################################################################
##
##
##
proc ModulesHelp { } {

   global software
   global version


   puts stderr "\tThis module provides support for ${software}"
   puts stderr "\n\tVersion: ${version}"

   puts stderr "\n\tSupport: h.french@garvan.org.au"
}

set     software        RNA-SeQC

set     version         1.1.7


module-whatis \
   "Setup environment for ${software} ${version} "

set env RNA-SeQC_HOME /share/ClusterShare/software/contrib/hugfre/RNA-SeQC/1.1.7

set prefix /share/ClusterShare/software/contrib/hugfre/RNA-SeQC/1.1.7

prepend-path    PATH            ${prefix}



