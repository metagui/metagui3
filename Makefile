###################################################################################
#                                                                                 #
#  Makefile for the compilation of metagui_util.f90                               #
#                                                                                 #
#  ------------ you may want to adapt these variables to your computer -----      #
#                                                                                 #
# ------------------------------------------------------------------------------- #
#                                                                                 #
#  This file is part of METAGUI.                                                  #
#                                                                                 #
#  METAGUI is free software: you can redistribute it and/or modify it under the   #
#  terms of the GNU General Public License as published by the Free Software      #
#  Foundation, either version 3 of the License, or (at your option) any later     #
#  version.                                                                       #
#                                                                                 #
#  METAGUI is distributed in the hope that it will be useful,  but WITHOUT ANY    # 
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR  #
#  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.    #
#                                                                                 #
#  You should have received a copy of the GNU General Public License along with   #
#  METAGUI.  If not, see <http://www.gnu.org/licenses/>.                          #
#                                                                                 #
###################################################################################

FC=gfortran
FFLAGS=-O3 -g
OMPFLAGS=-fopenmp
LDFLAGS=

SRC=src/metagui_util.f90 src/mk_microstates.f90 src/wham_bemeta.f90	\
    src/smart_clustering.f90 src/blas_dep.f src/dsyevx.f src/kinetic_2.f90


default: metagui_util.x


all: metagui_util.x metagui_util.exe


metagui_util.x: $(SRC)
	$(FC) $(OMPFLAGS) $(FFLAGS) -o $@ $^


metagui_util.exe: $(SRC)
	i686-w64-mingw32-gfortran $(OMPFLAGS)  $(FFLAGS) -s -static -o $@ $^


clean: 
	rm -f metagui_util.x  metagui_util.exe
