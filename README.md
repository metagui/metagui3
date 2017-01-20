Metagui 3
=========

METAGUI extends VMD with a graphical user interface that allows
constructing a thermodynamic and kinetic model of a given process
simulated by long-scale molecular dynamics. METAGUI can analyze both
unbiased and metadynamics-based simulations. The structures generated
during such a simulation are clustered together into a set of
microstates (i.e. structures with similar values of the collective
variables) and their relative free energies are then computed by a
weighted-histogram procedure.

The most relevant free energy wells are identified by diagonalization
of the rate matrix followed by a "smart clustering" analysis in
[Rodriguez-Laio 2014].  This procedure leads to a convenient
representation of the metastable states and long-time kinetics of the
system which can be compared with experimental data.

The tool allows to seamlessly switch between a collective variables
space representation of microstates and their atomic structure
representation, which greatly facilitates the set-up and analysis of
molecular dynamics trajectories.

METAGUI can read the output of the [PLUMED] 1.3 and 2.x engines, making
it compatible with a number of different molecular dynamics packages
like AMBER, NAMD, GROMACS and several others. 

[Rodriguez-Laio 2014]: http://www.sciencemag.org/content/344/6191/1492.abstract
[PLUMED]: http://www.plumed.org/


How to get it
-------------

Please find the latest release on GitHub at <https://github.com/metagui/metagui3>.



Installation instructions
-------------------------

First, the tool should be installed so that its root directory is
somewhere in TCL's "auto_path". Therefore, either move it into VMD's
plugin directory, or add the following instructions to your `.vmdrc`
file

        lappend auto_path /WHERE/YOU/EXTRACTED/THE/TOOL
        vmd_install_extension metagui3 metagui3::metagui "Analysis/Metagui 3.0"
        menu main on

Alternatively, you can set TCL's auto path via the `TCLLIBPATH`
environment variable.

You will also need to compile the FORTRAN code provided in the `src`
subdirectory (GNU gfortran is recommended): issue

        make
		
in the top directory. This should generate the executable
`metagui_util.x`. The plugin should automatically detect the location
of the executable. If it doesn't please copy it somewhere in your
PATH, or set the location in the "Advanced settings" tab.



Quickstart
----------

Just source the code as follows

        vmd -e run_metagui3.tcl

or, if you want to run from the DEMO_ALA3 directory

        TCLLIBPATH=.. vmd -e ../run_metagui3.tcl

See the read me file in that directory.



Authors
-------

The METAGUI 3.0 code builds upon
[METAGUI 2.0](http://dx.doi.org/10.1016/j.cpc.2011.08.020), written by
Xevi Biarnés, Fabio Pietrucci, Fabrizio Marinelli, Alessandro Laio.
Additions were made by Alex Rodriguez, Alessandro Laio and Toni
Giorgino.


License
-------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Additional software
-------------------

* Includes Paul Walton's scrollable frame code `sframe.tcl`, available
at http://wiki.tcl.tk/9223 .
* Includes the Tablelist 5.15 package, available from
http://www.nemethi.de/ .


Citation
--------

If you publish results obtained with METAGUI 3, please cite the
following paper

> Giorgino T, Laio A, Rodriguez A. METAGUI 3: a graphical user
> interface for choosing the collective variables in molecular
> dynamics simulations. (Submitted.)


The "smart clustering" algorithm is described in

> Rodriguez A, Laio A. Clustering by fast search and find of density
> peaks. Science. 2014;344(6191):1492-6. <http://www.sciencemag.org/content/344/6191/1492.abstract>

