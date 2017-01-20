# ------------------------------------------------------------------------------- #
# -------------------------------  M E T A G U I  ------------------------------- #
# --- A VMD Interface for Analyzing Long-Scale Molecular Dynamics Simulations --- #
# ------------------------------------------------------------------------------- #
#                                                                                 #
#  Copyright 2011-2016                                                            #
#										  #
# ------------------------------------------------------------------------------- #
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
# ------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------- #
# ---------------- G R A P H I C A L   U S E R   I N T E R F A C E -------------- #
# ------------------------------------------------------------------------------- #

package provide metagui3 3.0
package require tablelist
package require sframe


namespace eval ::metagui3:: {
    namespace export metagui

    # VARIOUS TEXTS AND CONSTANTS (Toni)
    # Message to show when all ok (constant)
    variable status_msg_ready "Ready."

    # CVs string for unbiased
    variable unloaded_string "(Still unloaded)"

    # Status message - automatically displayed
    variable status_msg $status_msg_ready

    # Some UTF-8 strings, required bc. windows assumes source in cp1252
    variable u8_dG "\u0394G";		   # \Delta G
    variable u8_rho "\u03c1";		   # \rho
    variable u8_delta "\u03b4";		   # \delta

    # Prefix for combined COLVAR files
    variable plumedgui_combined_prefix COMBINED

    # URLs - to revise
    variable citation_url "http://www.sciencedirect.com/science/article/pii/S0010465511003079"
    variable help_url "http://people.sissa.it/~laio/Publications/abstracts/65.php"

    # About message
    variable about_message "Current version: [package present metagui3]

Please cite the following paper below when publishing results
obtained with this tool:

FIXME"

    # Help for the decision plane procedure
    variable smart_clustering_help "Identify a region of isolated points in the high-$u8_rho, high-$u8_delta region of the decision plane.

Left mouse click on a point in the decision plane sets the cut-off for $u8_rho, right-click (or ctrl-left-click) sets it for $u8_delta.

See Rodriguez, A. & Laio, A. Clustering by fast search and find of density peaks. Science 344, 1492-1496 (2014)."

    # Help for the various tabs
    variable trajectory_tab_help "Set the structure and trajectories\
 to analyze. For each trajectory, a COLVAR file should be also\
 given.\nFor trajectories biased by metadynamics, set the\
 Biased property and provide the corresponding HILLS file."

    variable analysis_tab_help "(1) Choose a clustering algorithm to\
 derive microstates (they can visualized in the next pane).\n(2) Use\
 the \"thermodynamics\" tab to compute the the free energy of each\
 microstate.\nFor unbiased simulations, $u8_dG is proportional to\
 populations, whereas WHAM is used for biased runs.\n(3) Optionally, use\
 \"Clustering\" to aggregate microstates into kinetic basins."

    variable visualization_tab_help "Select 2 or 3 CVs to build projections, then select microstates or basins from the lists.
\"Show FES\" uses the states' $u8_dG for coloring or a third coordinate.
Keys F, D and E switch the 3D display mode; key P picks the microstate.   (Only available after analysis.)"

    # Help for Plumed-GUI integration
    variable plumedgui_button_help "You can define collective variables on-the-fly entering a PLUMED script in Plumed-GUI.\n
Selecting \"evaluate\" will run PLUMED's \"driver\" on each of the trajectory files listed.\n
The newly defined collective variables will be available for analysis (listed in green color).\n
Use the \"Merge\" function to merge the existing and new CVs into new COLVAR files."


    # Browse icon
    variable browse_icon [image create photo -file [file join $METAGUI_path icons openFolder.gif]]
    variable browse_icon_disabled [image create photo -file [file join $METAGUI_path icons openFolder_gray.gif]]


    # add VMD key short-cuts to control the display
    # show microstates (spheres)
    user add key f {metagui3::show_fes}
    # show structures
    user add key d {metagui3::show_system}
    # show both
    user add key e {metagui3::show_both}


    # TONI update everything
    proc update_lists {} {
	variable in1
	variable out1
	variable otp1
	list_inputs $in1.tr.files
	list_cvs    $in1.cv.list
	list_hills_conv .metagui3.hlf.nb.analys.nb.wham.conv.list
    }


    # 
    # TONI updates in the input table with the data in memory.
    # Replaces the old functions list_trajs, list_colvs, list_hills
    proc list_inputs {tbl} {
	variable nTRAJECTORIES
	variable TRAJ_INFO
	variable browse_icon
	variable browse_icon_disabled
	
	$tbl delete 0 end

	for {set i 0} {$i<$nTRAJECTORIES} {incr i} {
	    #	    set kkk [info exists TRAJ_INFO($i,bias)]
	    # if {$kkk==0} {set TRAJ_INFO($i,bias) 0}
	    set cvs $::metagui3::unloaded_string
	    if { [info exists TRAJ_INFO($i,hills_cvs)] && 
		 $TRAJ_INFO($i,hills_cvs)!="" }  {
		set cvs $TRAJ_INFO($i,hills_cvs)
	    }
    
	    $tbl insert end [list $TRAJ_INFO($i,traj_file) \
				 $TRAJ_INFO($i,temp) \
				 $TRAJ_INFO($i,colvar_file) \
				 $TRAJ_INFO($i,bias) \
				 $cvs \
				 $TRAJ_INFO($i,hills_file)]
	    		    
	    $tbl cellconfigure end,traj_file -image $browse_icon
	    $tbl cellconfigure end,file -image $browse_icon
	    $tbl cellconfigure end,hills -image $browse_icon

	    if {$TRAJ_INFO($i,bias) == 0} {
		$tbl cellconfigure end,biascv -editable no -foreground gray50 -selectforeground gray75
		$tbl cellconfigure end,hills -editable no -foreground gray50 -selectforeground gray75
		$tbl cellconfigure end,hills -image $browse_icon_disabled
	    } 

	    update

	    # If I knew how to color re-enabled rows, I could do this
#	    list_inputs_disable_unbiased $tbl
	}
    }

    # Edit done: copy table state into variables. May be optimized
    # avoiding the call to list_inputs (as long as hills is disabled)
    proc list_inputs_cellupdated {tbl} {
	variable nTRAJECTORIES
	variable TRAJ_INFO

	# A slight overkill - we update all the variables even though
	# we updated one cell only.
	set alldata [$tbl get 0 end]
	if {$nTRAJECTORIES != [llength $alldata]} {
	    error "BUG! nTRAJECTORIES=$nTRAJECTORIES while table has [llength $alldata] rows"
	}

	set i 0
	foreach row $alldata {
	    lassign $row TRAJ_INFO($i,traj_file) \
				 TRAJ_INFO($i,temp) \
				 TRAJ_INFO($i,colvar_file) \
				 TRAJ_INFO($i,bias) \
				 junk \
				 TRAJ_INFO($i,hills_file)
	    incr i
	}

	# TONI this is ugly- upon each edit, we empty and re-fill the
	# table list. It's done to make it easy to re-enable stuff - I
	# don't know how to recover the "non-disabled" background
	# color.
	list_inputs $tbl
    }

    # Add one row to the list of input files
    proc list_inputs_add {tbl} {
	variable nTRAJECTORIES
	variable TRAJ_INFO

	set n $nTRAJECTORIES
	incr nTRAJECTORIES
	set metagui3::TRAJ_INFO($n,traj_file) "trajectory.xtc"
	set metagui3::TRAJ_INFO($n,temp) 300
	set metagui3::TRAJ_INFO($n,colvar_file) "COLVAR"
	set metagui3::TRAJ_INFO($n,bias) 0
	set metagui3::TRAJ_INFO($n,biascv)	        "---"
	set metagui3::TRAJ_INFO($n,hills_file) "(none)"
	if {$n>0} {
	    # Try to provide sensible defaults copying from the previous trajectory
	    set p [expr {$n-1}]
	    set metagui3::TRAJ_INFO($n,traj_file)	$metagui3::TRAJ_INFO($p,traj_file) 
	    set metagui3::TRAJ_INFO($n,temp)		$metagui3::TRAJ_INFO($p,temp) 
	    set metagui3::TRAJ_INFO($n,colvar_file)	$metagui3::TRAJ_INFO($p,colvar_file) 
	    set metagui3::TRAJ_INFO($n,bias)		$metagui3::TRAJ_INFO($p,bias) 
	    set metagui3::TRAJ_INFO($n,hills_file)	$metagui3::TRAJ_INFO($p,hills_file)
	} 
	list_inputs $tbl
    }

    # Remove one row from the list of input files
    proc list_inputs_remove {tbl} {
	variable nTRAJECTORIES
	if {$nTRAJECTORIES>0} {
	    incr nTRAJECTORIES -1
	    list_inputs $tbl
	}
    }

    # Called when clicking on one of the folder icons in the input file tablelist
    proc list_inputs_browse {w x y} {
	variable TRAJ_INFO
	    # Check whether the mouse click occurred within a cell
	foreach {tbl x y} [tablelist::convEventFields $w $x $y] {}
	set cellIdx [$tbl containingcell $x $y]
	scan $cellIdx "%d,%d" row col
	if {$row < 0 || $col < 0} {
	    return ""
	}

	# Take an action if the clicked widget is the
	# label containing the image embedded in this cell
	if {$w eq [$tbl imagelabelpath $cellIdx]} {
	    set colname [$tbl columncget $col -name]
	    # Note: these tricky ifs use the && shortcut operator, so
	    # if nothing selected, no assignment is done
	    if {$colname=="traj_file" && [set nf [tk_getOpenFile]]!=""} {
		set TRAJ_INFO($row,traj_file) $nf 
	    } elseif {$colname=="file" && [set nf [tk_getOpenFile]]!=""} {
		set TRAJ_INFO($row,colvar_file) $nf 
	    } elseif {$colname=="hills" && $TRAJ_INFO($row,bias)==1 && [set nf [tk_getOpenFile]]!=""} { 
		set TRAJ_INFO($row,hills_file) $nf
	    }
	}
	update_lists

    }


    # proc list_inputs_disable_unbiased {tbl} {
    # 	variable nTRAJECTORIES
    # 	variable TRAJ_INFO
    # 	for {set i 0} {$i<$nTRAJECTORIES} {incr i} {
    # 	    if {$TRAJ_INFO($i,bias)==0} {
    # 		$tbl cellconfigure $i,biascv -foreground gray50 -selectforeground gray75 -editable no
    # 		$tbl cellconfigure $i,hills -foreground gray50 -selectforeground gray75 -editable no
    # 	    }
    # 	}
    # }


    # Update the CV tablelist 
    proc list_cvs {tbl} {
	variable nVARIABLES
	variable nRUNTIME
	variable CV_INFO_GRID
	variable COLVAR

#	parray ::metagui3::CV_INFO_GRID
	
	$tbl delete 0 end
	for {set i 1} {$i<=$nVARIABLES} {incr i} {
	    if { ![info exists CV_INFO_GRID($i,min)] || $CV_INFO_GRID($i,min)==""} \
		{ set CV_INFO_GRID($i,min) $COLVAR($i,min) }
	    if { ![info exists CV_INFO_GRID($i,max)] || $CV_INFO_GRID($i,max)==""} \
		{ set CV_INFO_GRID($i,max) $COLVAR($i,max) }
	    if { ![info exists CV_INFO_GRID($i,grid)] || $CV_INFO_GRID($i,grid)==""} \
		{ set CV_INFO_GRID($i,grid) 15 }
	    if { ![info exists CV_INFO_GRID($i,periodic)] || $CV_INFO_GRID($i,periodic)==""} \
		{ set CV_INFO_GRID($i,periodic) 0 }
	    if { ![info exists CV_INFO_GRID($i,use)] || $CV_INFO_GRID($i,use)==""} \
		{ set CV_INFO_GRID($i,use) 1 }

	    $tbl insert end [list "CV$i" \
				 $CV_INFO_GRID($i,type) \
				 $CV_INFO_GRID($i,min) \
				 $CV_INFO_GRID($i,max) \
				 $CV_INFO_GRID($i,grid) \
				 $CV_INFO_GRID($i,periodic) \
				 $CV_INFO_GRID($i,use) ]

	    # Visually mark runtime variables
	    set isRuntime [expr {$i>$nVARIABLES-$nRUNTIME}]
	    if {$isRuntime} {
		if { [expr {$i%2}] == 0 } {
		    set col PaleGreen
		} else {
		    set col LightGreen
		}
		$tbl rowconfigure end -background $col
	    }
	}

#	parray ::metagui3::CV_INFO_GRID

    }


    
    # React to interactive changes in the CV tablelist. Copy edits into CV_INFO_GRID.
    proc list_cvs_cellupdated {tbl} {
	variable nVARIABLES
	variable CV_INFO_GRID
	variable COLVAR

	set alldata [$tbl get 0 end]
	if {$nVARIABLES != [llength $alldata]} {
	    error "BUG! nVARIABLES=$nVARIABLES while table has [llength $alldata] rows"
	}

	# INFO GRID goes 1..N
	set i 1
	foreach row $alldata {
	    lassign $row junk \
		CV_INFO_GRID($i,type) \
		CV_INFO_GRID($i,min) \
		CV_INFO_GRID($i,max) \
		CV_INFO_GRID($i,grid) \
		CV_INFO_GRID($i,periodic) \
		CV_INFO_GRID($i,use) 
	    incr i
	}

        compute_active

    }



    # Updates the lower pane in Visualization/Microstates
    proc list_cvs_vis {cont} {
	variable nVARIABLES
	variable CV_INFO_GRID
	variable COLVAR

	catch {destroy $cont.me}
	sframe new $cont.me
	set myframe [sframe content $cont.me]

        frame $myframe.header
	label $myframe.header.id   -width 10 -text "CV #"
	label $myframe.header.name -width 12 -text "CV Name"
	label $myframe.header.plot           -text "Plot"
	pack $myframe.header.id $myframe.header.name $myframe.header.plot -side left -pady 2 
	pack $myframe.header -side top

	for {set i 1} {$i<=$nVARIABLES} {incr i} {
	    if { $CV_INFO_GRID($i,use) == 1 } {
		frame $myframe.$i 
		label $myframe.$i.id   -width 10 -text "CV$i: "
		entry $myframe.$i.type -width 12 -textvariable metagui3::CV_INFO_GRID($i,type) -background white
		checkbutton $myframe.$i.plot -variable metagui3::CV_INFO_GRID($i,useplot) \
		    -offvalue No -onvalue Yes -command metagui3::compute_active_plot
		pack $myframe.$i.id $myframe.$i.type $myframe.$i.plot -side left -pady 2 
		pack $myframe.$i -side top
	    }
	}

        sframe resize $cont.me
	pack $cont.me -side top
    }


    
    # Toni: not sure what this does
    proc list_hills_conv {myframe} {
	variable CV_INFO_GRID
	variable nTRAJECTORIES
	variable TRAJ_INFO

	puts "In list_hills_conv $myframe"
	catch {destroy $myframe.me}

	frame $myframe.me
	set NOTHING "---"

	set i 0
	for {set k 0} {$k<$nTRAJECTORIES} {incr k} {
	    if {$TRAJ_INFO($k,bias)==1} {
		### 
		#
		# A parlare, relazione tra CV_INFO; COLVARS & HILLS...
		#
		### 
		set j [expr " $k+1 "]
		if { ! [info exists CV_INFO_GRID($j,use)] } {
		    set CV_INFO_GRID($j,use) 0
		}
		if { $CV_INFO_GRID($j,use) == 1 } {
		    frame $myframe.me.$i 
		    entry $myframe.me.$i.hile -width 40 -background white \
			-textvariable metagui3::TRAJ_INFO([hills2traj $i],hills_file)
		    button $myframe.me.$i.showb -text "Plot" -height 1 \
			-command "metagui3::show_wham_graphs $i WHAM_RUN" 
		    pack $myframe.me.$i.hile  -side left  -expand yes -fill x
		    pack $myframe.me.$i.showb -side left  -expand no
		    pack $myframe.me.$i -side top -fill x
		}
		incr i
	    }
	}

	pack $myframe.me -side top
    }




    #
    # sets the default value of some variables. There are some more at
    # the beginning of metagui3.tcl (this function is probably redundant)

    proc bemeta_set_defaults {} {
	variable nTRAJECTORIES
	variable nVARIABLES
	
	variable nCV_ACTIVE
	variable CV_INFO_GRID
	variable CV_ACTIVE
	variable TRAJ_SKIP
	variable TRAJ_SKIP_BAFTI
	variable FES_MIN
	variable FES_MAX
	variable METAGUI_exe
	variable METAGUI_path
	variable DELTA
	variable GCORR
	variable TR_N_EXP
	variable REFES
	variable salign
	variable USE_ALL_WALKERS
	variable input_external
	variable cut_off
	variable sievlog
	variable sieving
	variable num_cluster
	variable cutting
	variable nconmin
	variable efdim
	variable FES_RANGE
	variable IGNORE_CONTROLS


	set IGNORE_CONTROLS 0
	set nTRAJECTORIES 0
	set nVARIABLES 0
	set nCV_ACTIVE 0

	set TRAJ_SKIP 1
	set TRAJ_SKIP_BAFTI 1

	set FES_MIN 0
	set FES_MAX 100
	set FES_RANGE 100

	set DELTA 4
	set GCORR 1
	set TR_N_EXP 5

	set cutting 1.0
	set nconmin 1
	
	set exe [expr {$::tcl_platform(platform)=="windows"?"exe":"x"}]
	set METAGUI_exe  [file join $METAGUI_path "metagui_util.$exe"]

	set REFES 0 

	set salign "noh"

	set USE_ALL_WALKERS 0

	set input_external GRID 
	set cut_off 1.0
	set sievlog 1
	set sieving 5000
	set num_cluster 100
	set efdim 1.0
    }

    # Hide/show the "sieve" command
    proc hidesieve {myframe} {
	variable sievlog
	variable sieving
	catch {destroy $myframe.me}
	frame $myframe.me
	switch $sievlog {
	    "0" {
		pack forget $myframe.me
	    }
	    "1" {
		frame $myframe.me.g
		label $myframe.me.g.label -text "Maximum distance matrix dimension: "
		entry $myframe.me.g.arg  -textvariable metagui3::sieving -background white -width 10
		pack $myframe.me.g.label -side left
		pack $myframe.me.g.arg   -side left -expand yes -fill x
		pack $myframe.me.g -side top -fill x
		pack $myframe.me -side top -fill x
	    }
	}
    }

    # hide arguments. Toni - these could be built once for all and
    # then shown/hidden by pack. Arguments are unnecessary.
    proc hideargs {myframe args_hide} {
	variable input_external
	variable cut_off
	variable sieving
	variable sievlog
	variable num_cluster
	variable out1
	catch {destroy $myframe.me}
	frame $myframe.me
	#    set sievlog "0"
	switch $args_hide {
	    "1" {
		set sievlog "0"
		pack forget $myframe.me
	    }
	    "2" {
		frame $myframe.me.g
		label $myframe.me.g.label -text "Cut off: "
		entry $myframe.me.g.arg  -textvariable metagui3::cut_off -background white -width 10
		pack $myframe.me.g.label -side left
		pack $myframe.me.g.arg   -side left -expand yes -fill x
		pack $myframe.me.g -side top -fill x
		frame $myframe.me.f
		checkbutton $myframe.me.f.arg -text "Sieving" -variable metagui3::sievlog \
                    -offvalue "0"  -onvalue "1" -command "metagui3::hidesieve $myframe.me.sieve"
		frame $myframe.me.sieve
		if { $sievlog == 1 } { hidesieve .metagui3.hlf.nb.analys.nb.cl.clust.ext.me.sieve }
		pack $myframe.me.f.arg -side left
		pack $myframe.me.f -side top -fill x
		pack $myframe.me.sieve -side top 
		pack $myframe.me -side top -fill x
	    }
	    "3" {
		frame $myframe.me.g
		label $myframe.me.g.label -text "Clusters: "
		entry $myframe.me.g.arg  -textvariable metagui3::num_cluster -background white -width 10
		pack $myframe.me.g.label -side left
		pack $myframe.me.g.arg   -side left -expand yes -fill x
		pack $myframe.me.g -side top -fill x
		frame $myframe.me.f
		checkbutton $myframe.me.f.arg -text "Sieving" -variable metagui3::sievlog \
                   -offvalue "0"  -onvalue "1" -command "metagui3::hidesieve $myframe.me.sieve"
		frame $myframe.me.sieve
		if { $sievlog == 1 } { hidesieve .metagui3.hlf.nb.analys.nb.cl.clust.ext.me.sieve }
		pack $myframe.me.f.arg -side left
		pack $myframe.me.f -side top -fill x
		pack $myframe.me.sieve -side top 
		pack $myframe.me -side top -fill x
	    }
	    "4" {		# TONI NEVER???
		frame $myframe.me.g
		label $myframe.me.g.label -text "cut off:"
		entry $myframe.me.g.arg  -textvariable metagui3::cut_off -background white -width 10
		pack $myframe.me.g.label -side left
		pack $myframe.me.g.arg   -side left -expand yes -fill x
		pack $myframe.me.g -side top -fill x
		frame $myframe.me.e
		label $myframe.me.e.label -text "efective dimension:"
		entry $myframe.me.e.arg  -textvariable metagui3::efdim -background white -width 10
		pack $myframe.me.e.label -side left
		pack $myframe.me.e.arg   -side left -expand yes -fill x
		pack $myframe.me.e -side top -fill x
		frame $myframe.me.f
		#    set sievlog "0"
		checkbutton $myframe.me.f.arg -text "sieving" -variable metagui3::sievlog \
		    -offvalue "0"  -onvalue "1" -command "metagui3::hidesieve $myframe.me.sieve"
		frame $myframe.me.sieve
		pack $myframe.me.f.arg -side left
		pack $myframe.me.f -side top -fill x
		pack $myframe.me.sieve -side top 
		pack $myframe.me -side top -fill x
	    }
	}
	#    if { $args_hide == 1 } {
	#          pack forget $myframe.me
	#       } else          {
	#          label $myframe.me.label -text "arguments:"
	#          entry  $myframe.me.arg  -textvariable metagui3::input_external  -background white -width 10
	#          pack $myframe.me.label -side left
	#          pack $myframe.me.arg   -side left -expand yes -fill x
	#          pack $myframe.me -side top -fill x
	#       } 
    }
    #
    #
    # hide traj frame
    proc hidetraj {myframe} {
	variable traj_hide
	catch {destroy $myframe.me}
	frame $myframe.me
	if { $traj_hide == 1 } {
	    pack forget $myframe.me
	    set traj_hide 0
	    $myframe.contain.htraj configure -text "(+)"
	} else          {
	    list_trajs $myframe.me
	    pack $myframe.me -side top -fill x
	    set traj_hide 1
	    $myframe.contain.htraj configure -text "(-)"
	} 
    }
    #
    # hide colvar frame
    proc hidecolvar {myframe} {
	variable colvar_hide
	catch {destroy $myframe.me}
	frame $myframe.me
	if { $colvar_hide == 1 } {
	    pack forget $myframe.me
	    set colvar_hide 0
	    $myframe.labels.hcolvar configure -text "(+)"
	} else          {
	    list_cvs $myframe.me
	    pack $myframe.me -side top
	    set colvar_hide 1
	    $myframe.labels.hcolvar configure -text "(-)"
	} 
    }
    #
    # hide hills frame
    proc hidehills {myframe} {
	variable hills_hide
	catch {destroy $myframe.me}
	frame $myframe.me
	if { $hills_hide == 1 } {
	    pack forget $myframe.me
	    set hills_hide 0
	    $myframe.labels.hhills configure -text "(+)"
	} else          {
	    list_hills $myframe.me
	    pack $myframe.me -side top -fill x
	    set hills_hide 1
	    $myframe.labels.hhills configure -text "(-)"
	} 
    }

    #
    # main GUI driver: creates the GUI
    # TONI TODO - use ttk as much as possible
    proc metagui {args} {
	variable KT
	variable T_CLUSTER
	variable T_FILL
	variable DELTA
	variable nVARIABLES
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable CV_INFO_GRID
	
	variable nTRAJECTORIES
	variable TRAJ_SKIP
	variable TRAJ_SKIP_BAFTI
	variable GRO_FILE
	
	variable TRAJ_INFO
	variable START_FOLDER  [pwd]
	variable GCORR
	variable TRAJ_molid
	variable TR_N_EXP 
	variable FES_molid
	# hide status
	variable traj_hide
	variable colvar_hide
	variable hills_hide
	variable salign
	variable FES_RANGE

	variable sievlog
	variable sieving

	variable CLUSTER_TYPE

	variable USE_ALL_WALKERS
	# variables for smart clustering
	variable RHO_CUT
	variable DELTA_CUT

	color Display {Background} white
	display projection   Orthographic
	display depthcue   off

	if { ![info exists nTRAJECTORIES] } {bemeta_set_defaults }

	variable w
	variable auto_show_clusters
	variable clusters_list
	variable auto_show_cluster_basin
	variable selected_cluster_basin_list

	# main window frame
	if { [winfo exists .metagui] } {
	    wm deiconify .metagui
	    return
	}
	set w [toplevel ".metagui3" -bg [ttk::style lookup . -background] -class Metagui2]
	wm title    $w "METAGUI [package present metagui3]"
	# allow .pippo_gui to expand with .
	grid columnconfigure $w 0 -weight 1
	grid rowconfigure $w 0 -weight 1
	wm geometry $w 825x700;	# instead it should expand properly
	# High Level Frame for tabs
	ttk::frame $w.hlf
	grid $w.hlf -column 0 -row 0 -sticky nsew
	# allow hlf to resize with window
	grid columnconfigure $w.hlf 0 -weight 1
	grid rowconfigure $w.hlf 0 -weight 1

	# notebook
	ttk::notebook $w.hlf.nb
	grid $w.hlf.nb -column 0 -row 0 -sticky nsew

	#---------------------------------------------------#
	#  Status line
	#---------------------------------------------------#
	ttk::frame $w.status
	grid .metagui3.status -sticky nswe 
	ttk::label $w.status.msg -textvariable metagui3::status_msg
	ttk::button $w.status.about -text "About" \
	    -command metagui3::help_about
	ttk::button $w.status.manual -text "Manual" \
	    -command "vmd_open_url $metagui3::help_url"
	ttk::button $w.status.citation -text "Citation" \
	    -command "vmd_open_url $metagui3::citation_url"
	pack $w.status.msg $w.status.about $w.status.manual $w.status.citation -side left
	pack $w.status.msg -expand 1

	#---------------------------------------------------#
	#  INPUT tab                            #
	#---------------------------------------------------#

	# build the frame, add it to the notebook
	ttk::frame $w.hlf.nb.simdat 
	$w.hlf.nb add $w.hlf.nb.simdat -text "Define inputs" -sticky nswe
	# allow frame to change width with window
	grid columnconfigure $w.hlf.nb.simdat 0 -weight 1 
	variable in1
	set in1 $w.hlf.nb.simdat

	labelframe $in1.menu 
	button $in1.menu.load -text "Load configuration" -command {metagui3::read_input [tk_getOpenFile]} -relief raised 
	button $in1.menu.save -text "Save configuration" -command {metagui3::write_input [tk_getSaveFile]} -relief raised 
	pack $in1.menu.load -side left
	pack $in1.menu.save -side left
	pack $in1.menu -side top

	# Trajectories/Colvars Section ----------------------------------

	labelframe $in1.tr -bd 2  -text "Trajectories" -padx 2m -pady 2m
	frame $in1.tr.gro
	label $in1.tr.gro.grol -text "Structure file:" -anchor e
	pack $in1.tr.gro.grol -side left
	button  $in1.tr.gro.grob -image $metagui3::browse_icon  -command {
	    if {[set tmp [tk_getOpenFile]]!=""} {set metagui3::GRO_FILE $tmp}
	} -pady 0 -padx 3
	pack $in1.tr.gro.grob -side left
	entry  $in1.tr.gro.grof -width 30 -textvariable metagui3::GRO_FILE  -background white
	pack $in1.tr.gro.grof -side left -expand yes -fill x
#	pack $in1.tr.gro  -side top -fill x
	grid $in1.tr.gro - -sticky ew

	# TONI The big tablelist for input. TG it would be nice if creation
	# and filling (list_inputs) could happen in the same
	# code. However, at this time TRAJ_INFO is empty.
	tablelist::tablelist $in1.tr.files \
	    -columntitles {
		"Trajectory file" 
		"Temperature" 
		"Colvar file"
		"Biased"
                "Bias CV"
		"Hills file"  } \
	    -stripebackground white  -instanttoggle 1  \
	    -yscrollcommand [list $in1.tr.vsb set]
	ttk::scrollbar $in1.tr.vsb -orient vertical -command [list $in1.tr.files yview]
	
	bind $in1.tr.files <<TablelistCellUpdated>> "metagui3::list_inputs_cellupdated %W"
	bind [$in1.tr.files bodytag] <Button-1> {::metagui3::list_inputs_browse %W %x %y}

	# These width and maxwidth do not seem to do anything
	$in1.tr.files columnconfigure 0 -editable true -name traj_file -align left -stretchable 1; # xtce
	$in1.tr.files columnconfigure 1 -editable true -name temp -align center ; # xtct
	$in1.tr.files columnconfigure 2 -editable true -name file -align left -stretchable 1; # colvare
	$in1.tr.files columnconfigure 3 -editable true -name bias -align center  -editwindow checkbutton  -formatcommand metagui3::bool2yesno; # colx
	$in1.tr.files columnconfigure 4 -editable false -name biascv -align center ; # NEW
	$in1.tr.files columnconfigure 5 -editable true -name hills -align left -stretchable 1; # hile
	#	pack $in1.tr.files $in1.tr.vsb -side left -fill both
	grid $in1.tr.files -sticky news
	grid $in1.tr.vsb -sticky ns -row 1 -column 1 
	grid columnconfigure $in1.tr 0 -weight 1
	# horizontal scrollbar in row 2 if desired



	frame $in1.tr.buttons
	button $in1.tr.buttons.addb -text "Add trajectory" -command "::metagui3::list_inputs_add $in1.tr.files"
	button $in1.tr.buttons.rmb -text "Remove trajectory" -command "::metagui3::list_inputs_remove $in1.tr.files"
	button $in1.tr.buttons.loadb -text "Load all"  -command "::metagui3::load_all"
	label  $in1.tr.buttons.loadinil -text " Start time:" -anchor e
	entry  $in1.tr.buttons.loadinie -width  6 -textvariable metagui3::T_CLUSTER -background white -validate key -validatecommand {
	    set ::metagui3::DELTA %P
	    return 1
	}
	label  $in1.tr.buttons.loadl -text " stride:" -anchor e
	entry  $in1.tr.buttons.loade -width  4 -textvariable metagui3::TRAJ_SKIP -background white

	pack $in1.tr.buttons.addb $in1.tr.buttons.rmb -side left
	pack $in1.tr.buttons.loade $in1.tr.buttons.loadl $in1.tr.buttons.loadinie $in1.tr.buttons.loadinil $in1.tr.buttons.loadb -side right
	grid $in1.tr.buttons - -sticky ew 

	pack $in1.tr -side top -fill both

	# end of Trajectoies/Colvars Section ----------------------------------


	labelframe $in1.plumedgui -text "Define new CVs on-the-fly with Plumed-GUI" -padx 2m -pady 2m
	pack [button $in1.plumedgui.help -bitmap questhead -command ::metagui3::plumedgui_help] -side right
	pack \
	    [button $in1.plumedgui.open -text "Open CV\neditor" -command "metagui3::plumedgui_open" ] \
	    [button $in1.plumedgui.eval_replace -text "Evaluate and\nreplace current" -command "metagui3::plumedgui_replace" ] \
	    [button $in1.plumedgui.eval_append -text "Evaluate and\nappend" -command "metagui3::plumedgui_append" ] \
	    [label $in1.plumedgui.spacer1 -text "      " ] \
	    [button $in1.plumedgui.rewrite -text "Merge existing\nand new CVs" -command "metagui3::plumedgui_combine" ] \
	    [label $in1.plumedgui.prefixlabel -text "into files named" ] \
	    [entry $in1.plumedgui.prefix -width 15 -background white -textvariable metagui3::plumedgui_combined_prefix ] \
	    [label $in1.plumedgui.spacer2 -text "      " ] \
	    -fill x -side left -expand 1
	pack $in1.plumedgui -side top -fill x

	# Collective Variables setup Section --------------------------------------

	labelframe $in1.cv -bd 2 -text "Collective variables (CVs)" -padx 2m -pady 2m

	tablelist::tablelist $in1.cv.list \
	    -columntitles { Id Name Min Max "Grid points" Periodic Use } \
	    -stretch all -stripebackground white -instanttoggle 1 \
	    -yscrollcommand [list $in1.cv.vsb set]
	ttk::scrollbar $in1.cv.vsb -orient vertical -command [list $in1.cv.list yview]
	
	bind $in1.cv.list <<TablelistCellUpdated>> "metagui3::list_cvs_cellupdated %W"

	$in1.cv.list columnconfigure 0 -editable false -name id -align center
	$in1.cv.list columnconfigure 1 -editable true -name name -align center
	$in1.cv.list columnconfigure 2 -editable true -name min -align center
	$in1.cv.list columnconfigure 3 -editable true -name max -align center
	$in1.cv.list columnconfigure 4 -editable true -name grid -align center
	$in1.cv.list columnconfigure 5 -editable true -name periodic -editwindow checkbutton \
	    -align center -formatcommand metagui3::bool2yesno
	$in1.cv.list columnconfigure 6 -editable true -name use -editwindow checkbutton \
	    -align center -formatcommand metagui3::bool2yesno
	pack $in1.cv.list  -side left -expand 1 -fill both
	pack $in1.cv.vsb -side left -fill y
	pack $in1.cv -side top -fill both

	# Help --------------
	ttk::label $in1.help -textvariable metagui3::trajectory_tab_help -padding 1m -justify center -anchor center 
	pack $in1.help -fill x -side bottom
	bind $in1.help <Configure>  { %W configure -wraplength [expr { %w - 20}] }


	# end of Collective Variables setup Section --------------------------------------

	# end of SETTINGS BOX --------------------------------------------------------------
	# ----------------------------------------------------------------------------------


	# ----------------------------------------------------------------------------------
	# ANALYSIS BOX ---------------------------------------------------------------------
	# build the frame, add it to the notebook
	ttk::frame $w.hlf.nb.analys 
	$w.hlf.nb add $w.hlf.nb.analys -text "Analyze"  -sticky nsew
	# allow frame to change width with window
	grid columnconfigure $w.hlf.nb.analys 0 -weight 1
	variable out1
	set out1 $w.hlf.nb.analys
	pack [ttk::notebook $out1.nb] -expand 1 -fill both

	# Cluster Analysis Section (.cl) ------------------------------------------------------
	# TAB MICROSTATES
	variable otp1
	set otp1 $w.hlf.nb.analys.nb
	$otp1 add [frame $otp1.cl] -text "Microstates"  -sticky nsew

	# Toni: this may get a face-lift. Fix TAG1 if so.
	frame $otp1.cl.clust 

	frame $otp1.cl.clust.comp -padx 10
	set statebutt disabled
	if {[info exists TRAJ_molid]} {set statebutt normal }
	button $otp1.cl.clust.comp.b -text "Find Microstates" \
                     -state $statebutt   -height 3 -command {
            metagui3::eval_with_status {Computing microstates...} {
		metagui3::do_clusters
		$metagui3::otp1.wham.ther.b configure -state normal 
		$metagui3::out2.cv.fes.show configure -state normal 
	    }
	}
	pack  $otp1.cl.clust.comp.b -side top -fill both -expand 1

	radiobutton $otp1.cl.clust.c -text "Grid" -variable metagui3::CLUSTER_TYPE -value grid \
            -command "set metagui3::input_external GRID ; metagui3::hideargs $otp1.cl.clust.ext 1"
	radiobutton $otp1.cl.clust.f -text "K-medoids (Kaufman and Rousseeuw)" -variable metagui3::CLUSTER_TYPE -value kmed \
	    -command "set metagui3::input_external KMEDOIDS ; metagui3::hideargs $otp1.cl.clust.ext 3"
	radiobutton $otp1.cl.clust.g -text "GROMOS (Daura et al.)" -variable metagui3::CLUSTER_TYPE -value grmc \
	    -command "set metagui3::input_external GROMACS ; metagui3::hideargs $otp1.cl.clust.ext 2"
	radiobutton $otp1.cl.clust.d -text "Read microstates" -variable metagui3::CLUSTER_TYPE -value read \
	    -command "metagui3::hideargs $otp1.cl.clust.ext 1"

	$otp1.cl.clust.c select

	frame $otp1.cl.clust.ext -padx 20

	pack $otp1.cl.clust.c -side top -anchor w
	pack $otp1.cl.clust.ext -side right -anchor nw -after $otp1.cl.clust.c
	pack $otp1.cl.clust.f -side top -anchor w
	pack $otp1.cl.clust.g -side top -anchor w
	pack $otp1.cl.clust.d -side top -anchor w
	pack $otp1.cl.clust.comp -side top -fill both -expand 1

	pack $otp1.cl.clust -side left -fill none

	#   pack $otp1.cl -side left  -fill none
	
	# end of Cluster Analysis Section (.cl) ---------------------------------------------------------


	# Wheighted-Histogram Analysis Section (.wham) ----------------------------------------------------
	$otp1 add [frame $otp1.wham] -text "Thermodynamics"  -sticky nsew
	labelframe $otp1.wham.ther -bd 2 -relief ridge -text "Compute $metagui3::u8_dG" -padx 4m -pady 4m
	frame $otp1.wham.ther.f1
	frame $otp1.wham.ther.f1.delta
	label $otp1.wham.ther.f1.delta.l -text "Delta" -anchor e
	entry $otp1.wham.ther.f1.delta.v -width 8 -textvariable metagui3::DELTA -background white
	pack $otp1.wham.ther.f1.delta.l -side left
	pack $otp1.wham.ther.f1.delta.v -side right
	pack $otp1.wham.ther.f1.delta -side top -fill x

	frame $otp1.wham.ther.f1.kt
	label $otp1.wham.ther.f1.kt.l -text "kT" -anchor e
	entry $otp1.wham.ther.f1.kt.v -width 8 -textvariable metagui3::KT -background white
	pack $otp1.wham.ther.f1.kt.l -side left
	pack $otp1.wham.ther.f1.kt.v -side right
	pack $otp1.wham.ther.f1.kt -side top -fill x

	frame $otp1.wham.ther.f1.tfill
	label $otp1.wham.ther.f1.tfill.l -text "Equilibration time" -anchor e
	entry $otp1.wham.ther.f1.tfill.v -width 8 -textvariable metagui3::T_FILL -background white
	pack $otp1.wham.ther.f1.tfill.l -side left
	pack $otp1.wham.ther.f1.tfill.v -side right
	pack $otp1.wham.ther.f1.tfill -side top -fill x
	pack $otp1.wham.ther.f1 -side top

	button $otp1.wham.ther.b -text "Run" -state disabled  -height 2 -pady 2 -command {
	    metagui3::eval_with_status {Computing free energies...} {
		metagui3::run_wham
#		$metagui3::otp1.smcl.cbas.run.c1 configure -state normal
		$metagui3::out2.cv.fes.b  configure -state normal
		$metagui3::otp1.smcl.scl.run.c1 configure -state normal
	    }
        }
	pack $otp1.wham.ther.b -side top -padx 10 -pady 1 -fill x -expand 1
	pack $otp1.wham.ther  -side left   -fill none -padx 4m



	# end of Wheighted-Histogram Analysis Section (.wham) ----------------------------------------------------

	labelframe $otp1.wham.conv -bd 2 -relief ridge -text "Check convergence" -padx 4m -pady 4m
	frame $otp1.wham.conv.list
	puts "list hills conv ini"
	list_hills_conv $otp1.wham.conv.list
	puts "list hills conv fin"
	pack $otp1.wham.conv.list -side top -fill x
	pack $otp1.wham.conv  -side left   -fill none -padx 4m
	#   pack $otp1.wham  -side left  -fill none -padx 4m


	$otp1 add [frame $otp1.smcl] -text "Clustering"  -sticky nsew

	# Kinetic Basins Section (.cbas) removed  ------------------------------------------------------------------

	# Smart clustering section
	labelframe $otp1.smcl.scl -bd 2 -relief ridge -text "Density-peaks Clustering" -padx 4m -pady 4m 

	frame  $otp1.smcl.scl.run; # TONI this frame is superfluous
	button $otp1.smcl.scl.run.c1 -text {Show decision plane} -state disabled   -height 3 -command {
                    metagui3::compute_active_plot 
                    $metagui3::otp1.smcl.scl.run.b configure -state normal 
                    metagui3::run_SAF
	}
	pack  $otp1.smcl.scl.run.c1 -side top  -pady 1 -fill x

	frame $otp1.smcl.scl.run.c2
	label $otp1.smcl.scl.run.c2.al -text "Cut off for $::metagui3::u8_rho :" -anchor e
	entry $otp1.smcl.scl.run.c2.at -width 15 -textvariable metagui3::RHO_CUT -background white
	pack $otp1.smcl.scl.run.c2.al -side left
	pack $otp1.smcl.scl.run.c2.at -side right
	pack $otp1.smcl.scl.run.c2 -side top -fill x

	frame $otp1.smcl.scl.run.c3
	label $otp1.smcl.scl.run.c3.al1 -text "Cut off for $::metagui3::u8_delta :" -anchor e
	entry $otp1.smcl.scl.run.c3.at1 -width 15 -textvariable metagui3::DELTA_CUT  -background white
	pack $otp1.smcl.scl.run.c3.al1 -side left
	pack $otp1.smcl.scl.run.c3.at1 -side right
	pack $otp1.smcl.scl.run.c3 -side top -fill x

	button $otp1.smcl.scl.run.b -text "Find kinetic basins via density-peaks clustering" -state disabled   -height 3 -command {
	    metagui3::eval_with_status {Building basins...} {
		metagui3::generate_input_smart
		$metagui3::otp1.smcl.scl.run.b configure -state disabled
		$metagui3::out2.cv.fes.bc configure -state normal 
	    }
	}
	pack $otp1.smcl.scl.run.b -side top   -fill x

	label $otp1.smcl.scl.run.h -pady 4m -text $::metagui3::smart_clustering_help
        bind $otp1.smcl.scl.run.h <Configure> { %W configure -wraplength [expr { %w - 4 }] }
	pack $otp1.smcl.scl.run.h -side top -fill x

	pack  $otp1.smcl.scl.run  -side top -fill both -expand 1
	pack $otp1.smcl.scl  -side left -fill both -expand 1 -padx 4m -pady 4m


	# Help --------------
	ttk::label $out1.help -textvariable metagui3::analysis_tab_help -padding 1m -justify center -anchor center 
	pack $out1.help -fill x -side bottom
	bind $out1.help <Configure>  { %W configure -wraplength [expr { %w - 20}] }


	# end of ANALYSIS BOX --------------------------------------------------------------
	# ----------------------------------------------------------------------------------


	# ----------------------------------------------------------------------------------
	# VISUALIZATION BOX ----------------------------------------------------------------
	# build the frame, add it to the notebook
	ttk::frame $w.hlf.nb.visual 
	$w.hlf.nb add $w.hlf.nb.visual -text "Visualize"  -sticky nsew
	# allow frame to change width with window
	grid columnconfigure $w.hlf.nb.visual 0 -weight 1
	variable out2
	set out2 $w.hlf.nb.visual

	# Clusters List Section (.cl_list) ------------------------------------------------------------------

	#    labelframe $out2.nb.cl_list -bd 2 -relief ridge -text "Microstates List" -padx 1m -pady 1m
	pack [ttk::notebook $out2.nb] -expand 1 -fill both
	$out2.nb add [frame $out2.nb.cl_list] -text "Microstates"  -sticky nsew
	frame  $out2.nb.cl_list.grid
	list_clusters
	pack $out2.nb.cl_list.grid -fill both -expand 1 -side top 

	frame $out2.nb.cl_list.checks
	checkbutton $out2.nb.cl_list.checks.autoshow -text "Show microstates in 3D view on select     " \
	    -variable metagui3::auto_show_clusters -offvalue 0 -onvalue 1 
	label $out2.nb.cl_list.checks.label -text "  Align to selection: " -anchor e
	entry $out2.nb.cl_list.checks.entry -textvariable metagui3::salign  -background white

	button $out2.nb.cl_list.checks.showb -text "Show" -height 1 -command {
	    metagui3::select_clusters "$metagui3::selected_clusters_list"
	    if {[info exists metagui3::FES_molid]} {
		metagui3::clean_fes
		metagui3::highlight_fes "$metagui3::selected_clusters_list"
	    } 
	}
	button $out2.nb.cl_list.checks.dumpb -text "Write PDB" -height 1 -command {
	    metagui3::write_clusters_frames_pdb "$metagui3::selected_clusters_list"
	} 
	pack $out2.nb.cl_list.checks.autoshow $out2.nb.cl_list.checks.label \
	     $out2.nb.cl_list.checks.entry $out2.nb.cl_list.checks.showb $out2.nb.cl_list.checks.dumpb -side left
	pack $out2.nb.cl_list.checks.entry -fill x -expand 1

	pack $out2.nb.cl_list.checks -side top -fill x 

	#    pack $out2.nb.cl_list -side top

	# end of Clusters List Section (.cl_list) ------------------------------------------------------------------

	# Basins List Section (.bs_list) ------------------------------------------------------------------

	#    labelframe $out2.nb.bs_list -bd 2 -relief ridge -text "Kinetic Basins List" -padx 1m -pady 1m
	$out2.nb add [frame $out2.nb.bs_list] -text "Basins"  -sticky nsew

	frame $out2.nb.bs_list.grid
	list_basins $out2.nb.bs_list.grid
	pack $out2.nb.bs_list.grid -side top -fill both -expand 1

	frame $out2.nb.bs_list.buttons
	checkbutton $out2.nb.bs_list.buttons.autoshow -text "Show basin center in 3D view on select     " \
	    -variable metagui3::auto_show_cluster_basin -offvalue 0 -onvalue 1 
	button $out2.nb.bs_list.buttons.showb -text "Show" -height 1 -command {metagui3::select_clusters "$metagui3::selected_cluster_basin_list"}
	pack $out2.nb.bs_list.buttons.autoshow -side left
	pack $out2.nb.bs_list.buttons.showb -side right
	pack $out2.nb.bs_list.buttons -side top -fill x

	#    pack $out2.nb.bs_list -side top -padx 4m

	# Collective Variables for visulaization --------------------------------------

	labelframe $out2.cv -bd 2 -relief ridge -text "FES and Basins plot" -padx 1m -pady 1m
	frame $out2.cv.list
	list_cvs_vis $out2.cv.list
	pack $out2.cv.list -side left
	
	labelframe $out2.cv.fes -text "View in 3D" -bd 2 -padx 4m -pady 4m
	button $out2.cv.fes.show -text "Microstates" -state disabled   -height 2 -command {
	    metagui3::compute_active_plot 
	    metagui3::show_basins $metagui3::var_list $metagui3::FES_RANGE MICROSTATES
	}
	pack  $out2.cv.fes.show -side top -fill x -expand 1

	frame $out2.cv.fes.al
	label $out2.cv.fes.al.b1 -text "$metagui3::u8_dG range: " -anchor e
	entry $out2.cv.fes.al.b2 -width 5 -textvariable metagui3::FES_RANGE -background white
	pack $out2.cv.fes.al.b1 -side left
	pack $out2.cv.fes.al.b2 -side left
	pack $out2.cv.fes.al -side top

	button $out2.cv.fes.b -text "Free Energy Surface" -state disabled  -height 2 -command {
	    metagui3::compute_active_plot 
	    metagui3::show_basins $metagui3::var_list $metagui3::FES_RANGE FES
	}
	pack $out2.cv.fes.b -side top -fill x -expand 1
	
	button $out2.cv.fes.bc -text "Basins" -state disabled  -height 2 -command {
	    metagui3::compute_active_plot 
	    metagui3::show_basins $metagui3::var_list $metagui3::FES_RANGE BASINS
	}
	pack $out2.cv.fes.bc -side top -fill x -expand 1

	pack $out2.cv.fes -side top -fill both -expand 1 
	pack $out2.cv -side top -fill both 


	# Help --------------
	ttk::label $out2.help -textvariable metagui3::visualization_tab_help -padding 1m -justify center -anchor center 
	pack $out2.help -fill x -side bottom
	bind $out2.help <Configure>  { %W configure -wraplength [expr { %w - 20}] }


	# end of VISUALIZATION BOX --------------------------------------------------------------
	# ----------------------------------------------------------------------------------

	# end of Clusters List Section (.cl_list) ------------------------------------------------------------------


	# ----------------------------------------------------------------------------------
	# ADVANCED ---------------------------------------------------------------------

	# build the frame, add it to the notebook
	ttk::frame $w.hlf.nb.adv 
	$w.hlf.nb add $w.hlf.nb.adv -text "Advanced settings"  -sticky nsew

	# allow frame to change width with window
	variable wadv $w.hlf.nb.adv
	grid columnconfigure $wadv 0 -weight 1
	grid columnconfigure $wadv 1 -weight 1

        grid [ ttk::label $wadv.ggl -text "Number of correlated frames in COLVARS (g in WHAM) " -anchor e ] \
 	     [ ttk::entry $wadv.ggv -width 10 -textvariable metagui3::GCORR ] 

        grid [ ttk::label $wadv.tnl -text "Minimum number of visits for a microstate " -anchor e ] \
             [ ttk::entry $wadv.tnv -width 10 -textvariable metagui3::TR_N_EXP ]

        grid [ ttk::label $wadv.rfl -text "Refine free energy " -anchor e ] \
	     [ ttk::checkbutton $wadv.rfv  -variable metagui3::REFES -offvalue 0  -onvalue 1 ] 

        grid [ ttk::label $wadv.kfl -text "Ignore convergence controls in FE calculation " -anchor e ] \
	     [ ttk::checkbutton $wadv.kfv  -variable metagui3::IGNORE_CONTROLS -offvalue 0  -onvalue 1 ]

        grid [ ttk::label $wadv.conl -text "Connectivity cutoff " -anchor e ] \
	     [ ttk::entry $wadv.conv -width 10 -textvariable metagui3::cutting ]

        grid [ ttk::label $wadv.alll -text "Use all frames (not recommended) " -anchor e ] \
	     [ ttk::checkbutton  $wadv.allv -variable metagui3::USE_ALL_WALKERS -offvalue 0 -onvalue 1  ]

        grid [ ttk::label $wadv.trnl -text "Minimum number of transitions to consider connected " -anchor e ] \
	     [ ttk::entry $wadv.trnv -width 10 -textvariable metagui3::nconmin  ]

        grid [ ttk::label $wadv.progl -text "Program location " -anchor e ] \
             [ ttk::frame $wadv.progf ]

        pack [ ttk::entry $wadv.progf.e -width 20 -textvariable metagui3::METAGUI_exe ] \
             [ ttk::button $wadv.progf.b -image $metagui3::browse_icon  -command {
 	                   if {[set tmp [tk_getOpenFile]]!=""} {set metagui3::METAGUI_exe $tmp}
              } ]  -side left

        grid [ ttk::labelframe $wadv.dd   -relief ridge -text "Diffusion matrix" ] -
	for { set i  1} {$i <= $nVARIABLES} { incr i } {
	    ttk::frame $wadv.dd.$i
	    for { set j  1} {$j <= $nVARIABLES} { incr j } {
		ttk::entry $wadv.dd.$i.$j -width 10 -textvariable metagui3::DD($i,$j)
		pack $wadv.dd.$i.$j -side left
	    }
	    pack $wadv.dd.$i -side top
	}

        foreach cell [grid slaves $wadv -column 0] {
             grid configure $cell -sticky e
        }
        foreach cell [grid slaves $wadv -column 1] {
             grid configure $cell -sticky w
        }

        grid configure $wadv.dd -sticky nsew

    }




    # CLUSTERS LIST ---------------------------------
    # scroll list data handler

    # creates the MICROSTATE section in the GUI. The frame passed as
    # an argument must exist.
    proc list_clusters {} {
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable clusters_tablelist
        variable out2

        set me $out2.nb.cl_list.grid

	catch {destroy $me.tbl $me.vsb $me.hsb}

	set cols {Id Size}
	for {set i 0} {$i < $nCV_ACTIVE} {incr i} {
	    lappend cols "CV$CV_ACTIVE($i)"
	}
	lappend cols "$metagui3::u8_dG" "Basin"

	set clusters_tablelist $me.tbl
	tablelist::tablelist $me.tbl -columntitles $cols \
	    -stretch all  -stripebackground white \
	    -selectmode extended  -labelcommand tablelist::sortByColumn  \
	    -xscrollcommand [list $me.hsb set] -yscrollcommand [list $me.vsb set] 
	ttk::scrollbar $me.vsb -orient vertical -command [list $me.tbl yview]
	ttk::scrollbar $me.hsb -orient horizontal -command [list $me.tbl xview]
	grid $me.tbl -sticky news
	grid $me.vsb -sticky ns -row 0 -column 1
	grid $me.hsb -sticky ew -row 1 -column 0
	grid columnconfigure $me 0 -weight 1
	grid rowconfigure $me 0 -weight 1
	
	for {set i 0} {$i < [$me.tbl columncount]} {incr i} {
	    $me.tbl columnconfigure $i -sortmode real -align right
	}

	# for {set i 0} {$i < $nCV_ACTIVE} {incr i} {
	#     $me.tbl columnconfigure [expr 2+$i] -name "cv$i"
	#     $me.tbl columnconfigure [expr 2+$i] -labelcommand tablelist::sortByColumn
	# }
	
	bind $me.tbl <<TablelistSelect>> "metagui3::clusters_multiple_sel %W"
    }


    #
    # fill in the clusters_tablelist
    proc clusters_datalist_fill {} {
	variable CLUSTERS
	variable CLUSTERS_ORDER
	variable NCL
	variable clusters_tablelist

	$clusters_tablelist delete 0 end

	for {set i 1} {$i <= $NCL} {incr i} {
	    set row {}
	    lappend row [format "%i" $i]; # id
	    lappend row [format "%i" $CLUSTERS($i,size)]; # size
	    set cv 0
	    foreach cc $CLUSTERS($i,center) {
		lappend row [format "%.3f" $cc]; # cv value at center
		incr cv
	    }
	    lappend row [format "%.1f" $CLUSTERS($i,fes)]
	    lappend row [format "%i" $CLUSTERS($i,basin)]
	    $clusters_tablelist insert end $row
	}
    }


    proc clusters_multiple_sel {widget} {
	variable clusters_tablelist
	variable auto_show_clusters
	variable selected_clusters_list
	variable FES_molid
	variable salign

	# convert rows to cluster ids (take the first column)
	set sel  [$widget  curselection]
	set selected_clusters_list {}
	foreach r $sel {
	    lappend selected_clusters_list [$widget getcells $r,0]
	}

	if {$auto_show_clusters==1} {
	    select_clusters $selected_clusters_list
	    if {[info exists FES_molid]} {
		clean_fes
		highlight_fes $selected_clusters_list
	    }
	}
	# puts "SELECTION $sel -> $selected_clusters_list"
    }





    # BASINS ---------------------------

    #
    # creates the BASINS section in the GUI
    proc list_basins {me} {
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable basins_tablelist

	catch {destroy $me.tbl $me.vsb $me.hsb}
	set basins_tablelist $me.tbl

	set cols [list Basin Size "Center microstate" $metagui3::u8_dG]

	tablelist::tablelist $me.tbl -columntitles $cols \
	    -stretch all -stripebackground white \
	    -selectmode extended  -labelcommand tablelist::sortByColumn  \
	    -xscrollcommand [list $me.hsb set] -yscrollcommand [list $me.vsb set] 
	ttk::scrollbar $me.vsb -orient vertical -command [list $me.tbl yview]
	ttk::scrollbar $me.hsb -orient horizontal -command [list $me.tbl xview]
	grid $me.tbl -sticky news
	grid $me.vsb -sticky ns -row 0 -column 1
	grid $me.hsb -sticky ew -row 1 -column 0
	grid columnconfigure $me 0 -weight 1
	grid rowconfigure $me 0 -weight 1

	for {set i 0} {$i < [$me.tbl columncount]} {incr i} {
	    $me.tbl columnconfigure $i -sortmode real -align right
	}
	
	bind $me.tbl <<TablelistSelect>> "metagui3::basins_multiple_sel %W"
    }


    #
    # defines the variables $basins_datalist
    proc basins_datalist_fill {} {
	variable CLUSTERS
	variable NCL
	variable KT
	variable NB
	variable BASINS
	variable START_FOLDER
	variable basins_tablelist

	$basins_tablelist delete 0 end

	for {set b 1} {$b <= $NB} {incr b} {
	    set row {}
	    if {$BASINS($b,size)>0} {
		lappend row  [format "%i" $b]
		lappend row  [format "%i" $BASINS($b,size)]
		lappend row  [format "%i" $BASINS($b,cluster_center)]
		lappend row  [format "%.1f" $BASINS($b,avrg_fes)]
	    }
	    $basins_tablelist insert end $row
	}

    }

    # scroll list data handler
    proc basins_multiple_sel {widget} {
	variable basins_datalist
	variable auto_show_cluster_basin
	variable selected_cluster_basin_list

	set sel [$widget curselection]
	set selected_cluster_basin_list ""

	foreach r $sel {
	    lappend selected_cluster_basin_list [$widget getcells $r,2]
	}

	if {$auto_show_cluster_basin==1} {
	    select_clusters "$selected_cluster_basin_list"
	}
    }




    # Toni: when is this called? Currently never. Perhaps should be called on any change to the CVs
    proc make_basinbutton_inactive {} {
	.metagui3.hlf.nb.analys.nb.wham.ther.b configure -state disabled
	# .metagui3.hlf.nb.analys.top.wham.b configure -state disabled
	# .metagui3.hlf.nb.analys.nb.smcl.cbas.run.c1 configure -state disabled
	# .metagui3.hlf.nb.analys.nb.smcl.cbas.run.b configure -state disabled
	list_cvs_vis .metagui3.hlf.nb.visual.cv.list
	list_hills_conv .metagui3.hlf.nb.analys.nb.wham.conv.list
    }





    ############################# Misc functions
    #
    # prints the "About" dialog
    proc help_about {args} {
	variable about_message

	set wabo .metagui_about
	catch {destroy $wabo}
	toplevel    $wabo
	wm title    $wabo "About METAGUI" 
	wm iconname $wabo "About" 
#	wm minsize  $wabo 220 200  

	label $wabo.l1 -font "bold" \
	    -text "METAGUI - A Unified Analysis Tool for Long-scale MD Simulations" 
	label $wabo.l2 -text $about_message
        button $wabo.ok -text "Got it" -command "destroy $wabo"

	pack $wabo.l1  $wabo.l2  $wabo.ok -side top -pady 10
    }



    #
    # changes the list of active biases (CV_ACTIVE_HILLS) when one cliks on a "use" button in the CV section
    proc compute_active_hills {} {
	variable TRAJ_INFO
	variable CV_ACTIVE_HILLS
	variable NCV_ACTIVE_HILLS

	set cv_list ""
	set NCV_ACTIVE_HILLS 0
	for { set h 0} {$h < [nHILLS] } { incr h } {
	    if { $TRAJ_INFO([hills2traj $h],hills_use) == 1 } {            
		foreach ia $TRAJ_INFO([hills2traj $h],hills_cvs) {
		    set good 1
		    foreach ib $cv_list {
			if { $ia == $ib } { set good 0 }
		    }
		    if { $good == 1 } { 
			lappend cv_list $ia 
			incr NCV_ACTIVE_HILLS
		    }
		}
	    }
	}
	set CV_ACTIVE_HILLS [lsort -integer $cv_list]
	puts $CV_ACTIVE_HILLS
    }


    #
    # changes the list of CVs used for plotting the free energy
    # (var_list) when one cliks on a "plot" button in the CV section
    proc compute_active_plot {} {
	variable CV_INFO_GRID
	variable nVARIABLES
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable var_list

	set ii 0
	set ntot 0
	set var_list ""
	for { set i  1} {$i <= $nVARIABLES} { incr i } {
	    # Toni: added "info exists" below because toggling the CV before analysis was giving error 
	    if {[info exists CV_INFO_GRID($i,useplot)] && $CV_INFO_GRID($i,useplot)=="Yes"} {
		incr ntot
		if {$CV_INFO_GRID($i,use)==0} {
		    set var_list "error"
		    puts "INFO: variable $i not active. It cannot be used for plotting"
                    return
		}
		lappend var_list $ii
		if {$ntot>3} {
		    puts "$var_list"
		    set var_list "error"
		    puts "INFO: can plot maximum 3-dimensional free energies"
                    return
		}
	    }
	    if {$CV_INFO_GRID($i,use)==1} {incr ii}
	}
        # TONI Add "F" and a "zero" CV up to 3 coordinates
        if {[llength $var_list]==1} {
            lappend var_list "F" "zero"
        } elseif {[llength $var_list]==2} {
	    lappend var_list "F"
	} else {
            puts "INFO: the number of active coordinates does not allow a plot (compute_active_plot)"
        }
    }



    # Misc functions
    # ------------------

    # Set status message, evaluate block, clear status even in case of error.
    proc eval_with_status {msg block} {
	variable status_msg
	variable status_msg_ready
	set status_msg $msg
	update
	if [ catch { 
	    uplevel $block
	} e ] {
	    set status_msg "Sorry, an error occurred. Please see console."
	    puts "$e"
	} else {
	    set status_msg $status_msg_ready
	}	    
    }



    # Boolean to Yes/No
    proc bool2yesno {v} {
	set r "??"
	switch $v {
	    0 {set r No}
	    1 {set r Yes}
	} 
	return $r
    }
	    

    # Return a temporary file name appropriate for the current OS
    # http://wiki.tcl.tk/772 . If we were in TCL 8.6 we could use
    # [file tempfile]
    proc tmpdir_func { } {
	global tcl_platform
	switch $tcl_platform(platform) {
	    unix {
		set tmpdir /tmp   ;# or even $::env(TMPDIR), at times.
	    } macintosh {
		set tmpdir $::env(TRASH_FOLDER)  ;# a better place?
	    } default {
		set tmpdir [pwd]
		catch {set tmpdir $::env(TMP)}
		catch {set tmpdir $::env(TEMP)}
	    }
	}
	return $tmpdir
    }


    # Load everything
    proc load_all {} {
	metagui3::eval_with_status {Loading COLVAR and HILLS files...}   {
	    metagui3::load_data
	}
	metagui3::eval_with_status {Loading trajectories...} {
	    metagui3::load_trajs
	}
	.metagui3.hlf.nb.analys.nb.cl.clust.comp.b configure -state normal; # TAG1
	metagui3::update_lists
    }


    
    #
    # PLUMED-GUI 
    # ---------------


    # Open "Plumed-GUI cv analysis tool". Needs a recent version.
    proc plumedgui_open {} {
	if [winfo exists .plumed] {
	    wm deiconify .plumed
	    raise .plumed
	} else {
	    ::Plumed::plumed
	}
    }



    # Evaluate and load in place of the old COLVARS. Set everything as unbiased
    proc plumedgui_replace {} {
        variable nTRAJECTORIES
        variable TRAJ_INFO
        plumedgui_prepare
        set rcv_basename "RUNTIME_COLVAR"
	plumedgui_eval $rcv_basename
	for {set i 0} {$i<$nTRAJECTORIES} {incr i} {
            set TRAJ_INFO($i,colvar_file)  [format "%s.%03d" $rcv_basename $i]
            set TRAJ_INFO($i,bias) 0
            set TRAJ_INFO($i,hills_cvs) $::metagui3::unloaded_string
            unset TRAJ_INFO($i,runtime_colvar)
	}
        update_lists
    }


    # Evaluate and load new CVs 
    proc plumedgui_append {} {
        plumedgui_prepare

	set rcv_basename [file join [tmpdir_func] "RUNTIME_COLVAR.[pid]"]
	plumedgui_eval $rcv_basename

	eval_with_status "Loading newly computed COLVAR files..." {
	    load_runtime
	}
    }


    # Sanity checks for both replace and append
    proc plumedgui_prepare {} {
	variable nTRAJECTORIES
	if {![winfo exists .plumed]} {
	    error "Please open the Plumed-GUI CV editor first"
	}
	if {$nTRAJECTORIES<1} {
	    error "Please load trajectories first"
	}
    }


    # Auxiliary function. Launch plumed on the script currently
    # displayed in the "Plumed-GUI cv analysis tool". Write results in
    # temporary files.
    proc plumedgui_eval {rcv_basename} {
	variable nTRAJECTORIES
	variable TRAJ_INFO
	variable GRO_FILE

	set cm [mol new $GRO_FILE waitfor all]
	# Iterate over trajectories
	for {set i 0} {$i<$nTRAJECTORIES} {incr i} {
	    set statusmsg [format "Running PLUMED... (trajectory %d/%d)" [expr $i+1] $nTRAJECTORIES]
	    eval_with_status $statusmsg {
		set rcv [format "%s.%03d" $rcv_basename $i]
		set TRAJ_INFO($i,runtime_colvar) $rcv
		animate delete all
		mol addfile "$TRAJ_INFO($i,traj_file)" waitfor all
		::Plumed::do_compute $rcv; # Requires a recent VMD + Plumed-GUI. May fail.
	    }
	}
	mol delete $cm
    }


    # Auxiliary function.  Parse newly-created RUNTIME_COLVAR
    # files. THIS WILL INCREASE nVARIABLES!
    proc load_runtime {} {
	variable TRAJ_INFO;	# r
	variable nTRAJECTORIES;	# r
	variable nVARIABLES;	# w
	variable nRUNTIME;	# w
	variable CV_INFO_GRID;	# w
	variable COLVAR;	# w

	variable in1;		# gui stuff, needs refactoring
	variable out2;		# gui stuff
	
	# Probe the first only; supposedly equal.
	set probe_colvar $TRAJ_INFO(0,runtime_colvar)
	set rvarlist [parse_colvars_metadata $probe_colvar]
	set nRUNTIME [llength $rvarlist]

	# Now we essentially do what proc load_data does, but on
	# colvars between nVARIABLES+1 nVARIABLES+nRUNTIME (both
	# included).
	set i [expr $nVARIABLES+1]
	foreach varname $rvarlist {
	    set CV_INFO_GRID($i,type) $varname
	    set COLVAR($i,min) 9999999
	    set COLVAR($i,max) -9999999
	    incr i
	}

	# read and load into memory each colvar file
	for { set c 0} {$c < $nTRAJECTORIES } { incr c } {
	    puts ""
	    puts "reading runtime COLVAR file $TRAJ_INFO($c,runtime_colvar)" 
	    puts ""
	    read_colvars $c 1;	# note the runtime flag
	}

	set nVARIABLES [expr {$nVARIABLES+$nRUNTIME}]
	
	list_cvs $in1.cv.list
	list_cvs_vis $out2.cv.list
	#XXX
	list_hills_conv .metagui3.hlf.nb.analys.nb.wham.conv.list
	list_clusters 
	# automatic build the default difusion matrix
	build_DD
	
    }


    # Show help
    proc plumedgui_help {} {
         variable plumedgui_button_help
#         tk_dialog .metagui3.plumedgui_help "Using Plumed-GUI" $plumedgui_button_help questhead 0 "Ok"
         tk_messageBox -message $plumedgui_button_help \
                -title "Help" -parent .metagui3 -icon info -type ok
    }




    # Combine existing and new (runtime) COLVAR files. Headers need some care.
    proc plumedgui_combine {} {
	variable plumedgui_combined_prefix
	variable TRAJ_INFO
	variable nTRAJECTORIES

	eval_with_status "Merging old and newly-computed COLVAR files..." {
	    for {set tr 0} {$tr<$nTRAJECTORIES} {incr tr} {
		set in1 $TRAJ_INFO($tr,colvar_file)
		set in2 $TRAJ_INFO($tr,runtime_colvar)
		set out [format "%s.%03d" $plumedgui_combined_prefix $tr]
		puts "Merging $in1 + $in2 --> $out"
		combine_colvars $in1 $in2 $out
	    }
	}
    }

    
    # Auxiliary function. Combine two COLVAR files. Assumes FIELDS
    # header, and ACTIVE header on the first one. Ignores ACTIVE and
    # the time column of the second. Puts vbias last, if present. CV
    # name conflicts are not handled.
    proc combine_colvars {in1f in2f outf} {
	set in1 [open $in1f r]
	set in2 [open $in2f r]
	set out [open $outf w]

	# Fields
	gets $in1 fi1
	gets $in2 fi2

	set biasPresent 0;	# Whether the 1st input has a bias column
	if { [lindex $fi1 end] == "vbias" } {
	    set biasPresent 1
	}
	if { [lindex $fi2 end] == "vbias" } {
	    error "Runtime colvars should not have a bias column!"
	}

	# Drop hash
	set fi1 [lreplace $fi1 0 0]
	
	# Drop hash, FIELDS, and time headers
	set fi2 [lreplace $fi2 0 2]
	puts -nonewline $out "#! "
	puts $out [combine_colvars_mix $fi1 $fi2 $biasPresent]

	# ACTIVE HEADER on 1st file, if present - copy it
	gets $in1 line1
	if [regexp {ACTIVE} $line1] {
	    puts $out $line1
	    gets $in1 line1
	} 
	    
	# ACTIVE HEADER on 2nd file, if present - drop it
	gets $in2 line2
	if [regexp {ACTIVE} $line2] {
	    gets $in2 line2
	} 

	# From now on we assume equal number of remaining
	# lines.
	while {![eof $in1] && ![eof $in2]} {
	    set line2 [lreplace $line2 0 0]; # remove time
	    puts $out [combine_colvars_mix $line1 $line2 $biasPresent]
	    gets $in1 line1
	    gets $in2 line2
	}

        # Check for equal lenghts. Boths eofs should be 1 now
        if { [eof $in1] != [eof $in2] } {
            puts "WARNING: Colvars files have different number of timesteps, and have been truncated to the shortest"
        }

	close $in1
	close $in2
	close $out
    }

    # Auxiliary function.  Append elements of l2, except the first, to
    # the end (beforelast=0) or before the last element (beforelast=1)
    # of l1
    proc combine_colvars_mix {l1 l2 {beforelast 1}} {
	if {$beforelast==0} {
	    set where end
	} else {
	    set where {end-1}
	}
	return [linsert $l1 $where {*}$l2]
    }



}
