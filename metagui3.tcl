package provide metagui3 3.0

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


package require plumed

namespace eval ::metagui3:: {
    namespace export metagui

    #
    # declaration of plugin-specific variables
    #

    # path to metagui_util.x (precompiled tasks)
    # TG windows will also need auto_execok, to be fixed - see proc bemeta_set_defaults
    variable METAGUI_exe 

    # path to this file (must be done here)
    variable METAGUI_path [file normalize [file dirname [info script]]]

    # path to current working directory
    variable START_FOLDER

    # CLUSTERS(*,*,*) matrix containing microstates centers, free energy, basins, ...
    variable CLUSTERS
    # BASINS(*,*,*) matrix containing basins attractors, free energy, ...
    variable BASINS
    # FES(*,*) matrix containing Free Energy of each micrstate
    variable FES
    # total number of microstates
    variable NCL
    # total number of basins
    variable NB
    # difusion matrix
    variable DD
    # boolean to refine free energy during kinetic analysis
    variable REFES

    # variables for internal use only
    variable CLUSTERS_FES_list
    variable CLUSTERS_bin
    variable CLUSTERS_center
    variable CLUSTERS_size
    variable clusters_LIST
    variable ncLIST
    variable current_cluster
    variable CLUSTERS_NCL
    variable CONNECTED_GROUPS
    variable NG
    variable FES_MIN
    variable FES_MAX
    variable reprFES_MAX
    variable nREPSinFES
    variable pickingON
    # variables for smart clustering
    variable RHO_CUT
    variable DELTA_CUT


    # maximum number of structures to represent for each microstate (to avoid memory outage in VMD display)
    variable MAX_FRAMES_PER_CLUSTER
    # hardcoded!!
    set MAX_FRAMES_PER_CLUSTER 100

    # VMD mol id containing the trajectory
    variable TRAJ_molid
    # VMD mol id containing the microstates representation (spheres)
    variable FES_molid

    # total number of active CV for microstates analysis
    variable nCV_ACTIVE
    # list of active CV for microstates analysis
    variable CV_ACTIVE
    # grid size for each CV for microstates analysis
    variable CV_INFO_GRID
    # COLVAR(*,*,*) matrix containing the values of each CV at each time for each COLVAR file
    variable COLVAR
    # total number of trajectories loaded
    variable nTRAJECTORIES
    # time after which to load the trajectories
    variable T_CLUSTER   1
    # time after wich to perform the WHAM analysis
    variable T_FILL      1
    # WHAM analysis parameters
    variable DELTA       4
    variable KT          2.4943
    variable GCORR
    variable TR_N_EXP 
    # skip TRAJ_SKIP frames when loading trajectories
    variable TRAJ_SKIP
    variable TRAJ_SKIP_BAFTI
    # structure file for loading trajectory XTC files
    variable GRO_FILE
    # information on each COLVAR file
    variable TRAJ_INFO

    # total number of active variables
    variable nVARIABLES
    variable var_list

    # total number of frames loaded from COLVAR files
    variable TOTAL_COLVAR_FRAMES
    # total number of frames loaded from trajectory XTC files
    variable TOTAL_COLVAR_STRUCTURES
    # total number of frames in trajectory XTC files
    variable TOTAL_TRAJ_FRAMES
    # total number of frames in VMD
    variable FRAMES

    # GRID dimensions for each CV
    variable DIMENSIONS

    # GRID space for each active CV
    variable CLUSTER_SPACE
    variable CLUSTER_SPACE_dim
    variable CLUSTER_SPACE_vol
    # list of active CV
    variable myVARS
    # list of CVs active in HILLS files
    variable CV_ACTIVE_HILLS
    # total number of CVs active in HILLS files
    variable NCV_ACTIVE_HILLS

    # ordered list of microstates
    variable CLUSTERS_ORDER

    # clustering method type (grid, external, read)
    variable CLUSTER_TYPE
    variable sievlog

    # use all colvar frames for clustering analysis, and not only those biased by active hills in use
    variable USE_ALL_WALKERS
    # ignore convergence controls in Free Energy calculation
    variable IGNORE_CONTROLS

    # plotting with Multiplot
    variable ploth
    variable winf .metagui_plot


    #
    # read input file and store information in variable variables
    proc read_input {{INPUT_file "metagui.in"}} {
	variable KT
	variable T_CLUSTER
	variable T_FILL
	variable GCORR
	variable TR_N_EXP
	variable DELTA
	variable nVARIABLES
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable CV_INFO_GRID
	
	variable nTRAJECTORIES
	variable TRAJ_SKIP
	variable TRAJ_SKIP_BAFTI
	variable GRO_FILE
	variable CLUSTER_TYPE
	variable sieving
	variable sievlog
	variable cut_off
	variable num_cluster
	variable FES_RANGE
	
	
	variable TRAJ_INFO
	variable FES_MIN
	variable FES_MAX
	variable METAGUI_exe
	variable START_FOLDER

	set TRAJ_SKIP 10
	set TRAJ_SKIP_BAFTI 1

	set h 0
	set c 0
	if {[file exists $INPUT_file]} {
	    set in [open $INPUT_file r]
	    while {[gets $in line] != -1} {
		if [regexp {^#} $line] {
		    continue;	# skip comments
		}
		set KEY [lindex $line 0]
		switch $KEY {
		    "METAGUI_EXE" {
			set METAGUI_exe [lindex $line 1]
		    }
		    "KT" {
			set KT [lindex $line 1]
		    }
		    "T_CLUSTER" {
			set T_CLUSTER [lindex $line 1]
		    }
		    "T_FILL" {
			set T_FILL [lindex $line 1]
		    }
		    "GCORR" {
			set GCORR [lindex $line 1]
		    }
		    "TR_N_EXP" {
			set TR_N_EXP [lindex $line 1]
		    }
		    "DELTA" {
			set DELTA [lindex $line 1]
		    }
		    # There's no need to provide this anymore. nVARIABLES is get from parse_colvars
		    #        "NCV" {
		    #          set nVARIABLES [lindex $line 1]
		    #        }
		    "ACTIVE" {
			set nCV_ACTIVE [lindex $line 1]
			for { set i  0} {$i < $nCV_ACTIVE } { incr i } {
			    set cv [lindex $line [expr "2 + $i"]]
			    set CV_ACTIVE($i) $cv
			    set CV_INFO_GRID($cv,use) 1
			}
		    }
		    "CVGRID" {
			set idum [lindex $line 1]
			set CV_INFO_GRID($idum,min) [lindex $line 2] 
			set CV_INFO_GRID($idum,max) [lindex $line 3] 
			set CV_INFO_GRID($idum,grid) [lindex $line 4] 
			set tmp [lindex $line 5]
			if {$tmp=="PERIODIC"} { 
			    set CV_INFO_GRID($idum,periodic) 1 
			    set CV_INFO_GRID($idum,type) " "
			} else {
			    set CV_INFO_GRID($idum,type) [lindex $line 5]
			}
			set tmp [lindex $line 6]
			set CV_INFO_GRID($idum,periodic) 0
			if {$tmp=="PERIODIC"} { set CV_INFO_GRID($idum,periodic) 1 }
			set kkk [info exists CV_INFO_GRID($idum,use)]
			if {$kkk==0} {set CV_INFO_GRID($idum,use) 0}
		    }
		    "HILLS_FILE" {
			puts "WARNING: HILLS_FILE is obsolete in input."
		    }
		    "COLVAR_FILE" {
			set TRAJ_INFO($c,colvar_file) [lindex $line 1] 
			set TRAJ_INFO($c,traj_file) [lindex $line 2] 
			set TRAJ_INFO($c,temp) [lindex $line 3]
			set TRAJ_INFO($c,bias) [lindex $line 4]
			if {$TRAJ_INFO($c,bias)==1} {
			    set TRAJ_INFO($c,hills_file) [lindex $line 5]
			} else {
			    set TRAJ_INFO($c,hills_file) ""
			}
			incr c
		    }
		    "TRAJ_SKIP" {
			set TRAJ_SKIP [lindex $line 1]
		    }
		    "TRAJ_SKIP_BAFTI" {
			set TRAJ_SKIP_BAFTI [lindex $line 1]
		    }
		    "GRO_FILE" {
			set GRO_FILE [lindex $line 1]
		    }
		    "CLUSTER_TYPE" {
			set CLUSTER_TYPE [lindex $line 1]
		    }
		    "SIEVING" {
			set sievlog [lindex $line 1]
		    }
		    "NUMBER_STRUCTURES_SIEVING" {
			set sieving [lindex $line 1]
		    }
		    "MEDOIDS_CLUSTERS" {
			set num_cluster [lindex $line 1]
		    }
		    "GROMACS_CUTOFF" {
			set cut_off [lindex $line 1]
		    }
		    "FES_RANGE" {
			set FES_RANGE [lindex $line 1]
		    }
		    default {
			puts "WARNING! mis understood line:"
			puts "  $line"
		    }
		}
	    }
	    #   compute_active_hills
	} else {
	    puts ""
	    puts "WARNING: input file $INPUT_file not found"
	    puts ""
	    return "input file $INPUT_file not found"
	}

	set nTRAJECTORIES $c

	puts ""
	puts "INPUT file ($INPUT_file) read succesfully"
	puts ""

	update_lists
	if {$CLUSTER_TYPE == "kmed"} {hideargs .metagui3.hlf.nb.analys.nb.cl.clust.ext 3}
	if {$CLUSTER_TYPE == "grmc"} {hideargs .metagui3.hlf.nb.analys.nb.cl.clust.ext 2}

    }

    #
    # write configuration into a file
    proc write_input {{NEWINPUT_file "cluster_new.in"}} {
	variable KT
	variable T_CLUSTER
	variable T_FILL
	variable GCORR
	variable TR_N_EXP
	variable DELTA
	variable nVARIABLES
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable CV_INFO_GRID
	
	variable nTRAJECTORIES
	variable TRAJ_SKIP
	variable TRAJ_SKIP_BAFTI
	variable GRO_FILE
	variable FES_RANGE
	variable CLUSTER_TYPE
	variable REPR
	variable NR
	
	
	variable TRAJ_INFO
	variable FES_MIN
	variable FES_MAX
	variable METAGUI_exe
	variable START_FOLDER
	variable sieving
	variable sievlog
	variable cut_off
	variable num_cluster

	set out [open $NEWINPUT_file w]

	puts $out "# METAGUI_EXE $METAGUI_exe"
	puts $out "KT $KT"
	# TONI FIXME
	# for { set h  0} {$h < [nHILLS] } { incr h } {
	#     set kkk [info exists TRAJ_INFO([hills2traj $h],hills_use)]
	#     if { $kkk == 0 } { set TRAJ_INFO([hills2traj $h],hills_use) 0 }
	#     puts $out "HILLS_FILE $TRAJ_INFO([hills2traj $h],hills_file) $TRAJ_INFO([hills2traj $h],hills_use) $HILLS_INFO($h,traj)"
	# }
	puts $out "CLUSTER_TYPE $CLUSTER_TYPE"
	puts $out "FES_RANGE $FES_RANGE"
	puts $out "GRO_FILE $GRO_FILE"
	for { set c 0} {$c < $nTRAJECTORIES } { incr c } {
	    puts $out [list COLVAR_FILE \
			   $TRAJ_INFO($c,colvar_file) \
			   $TRAJ_INFO($c,traj_file) \
			   $TRAJ_INFO($c,temp) \
			   $TRAJ_INFO($c,bias) \
			   $TRAJ_INFO($c,hills_file)]
	}
	puts $out "TRAJ_SKIP $TRAJ_SKIP"
	for {set i 1} {$i < [expr "$nVARIABLES+1"]} {incr i} {
	    puts -nonewline $out  "CVGRID  $i $CV_INFO_GRID($i,min) $CV_INFO_GRID($i,max) $CV_INFO_GRID($i,grid) $CV_INFO_GRID($i,type)"
	    if { $CV_INFO_GRID($i,periodic)==1 } { puts  -nonewline $out " PERIODIC" }
	    puts $out ""
	}
	puts -nonewline $out  "ACTIVE $nCV_ACTIVE  "
	for { set i  0} {$i < $nCV_ACTIVE } { incr i } {
	    puts -nonewline $out  "$CV_ACTIVE($i) "
	}
	puts $out ""
	puts $out "T_CLUSTER $T_CLUSTER"
	puts $out "T_FILL $T_FILL"
	puts $out "DELTA $DELTA"
	puts $out "GCORR $GCORR"
	puts $out "TR_N_EXP $TR_N_EXP"
	puts $out "SIEVING $sievlog"
	puts $out "NUMBER_STRUCTURES_SIEVING $sieving"
	puts $out "MEDOIDS_CLUSTERS $num_cluster"
	puts $out "GROMACS_CUTOFF $cut_off"
	close $out
    }


    

    #
    # load trajectories from XTC files
    proc load_trajs {} {
	variable nTRAJECTORIES
	variable COLVAR
	variable TRAJ_SKIP
	variable GRO_FILE
	variable TOTAL_TRAJ_FRAMES
	variable TRAJ_molid
	variable TOTAL_COLVAR_STRUCTURES
	variable START_FOLDER
	cd $START_FOLDER
	variable TRAJ_INFO

	set prev_TOTAL_TRAJ_FRAMES 0

	# load structure file first and delete the coordinates
	mol new $GRO_FILE first 0 last -1 step 1  filebonds 1 autobonds 1 waitfor all
	set TRAJ_molid [molinfo top get id]
	animate delete beg 0 end 0 $TRAJ_molid

	for { set c 0} {$c < $nTRAJECTORIES } { incr c } {
	    set colvar_file $TRAJ_INFO($c,colvar_file)
	    set colvar_traj_file $TRAJ_INFO($c,traj_file)
	    # check wether COLVAR file has been loaded first
	    if {! [info exists COLVAR($c,traj_frames)]} {
		error "ERROR! trajectories incorrectly loaded because colvar data was not loaded previously. Please, use 'load_data' first."
	    }
	    puts ""
	    puts "loading TRAJECTORY file $colvar_traj_file"
	    puts "starting from frame $COLVAR($c,first_traj_frame_to_be_loaded) and skiping every $TRAJ_SKIP frames"
	    mol addfile $colvar_traj_file first $COLVAR($c,first_traj_frame_to_be_loaded) last -1 step $TRAJ_SKIP filebonds 1 autobonds 1 waitfor all
	    set TOTAL_TRAJ_FRAMES [molinfo $TRAJ_molid get numframes]
	    set nFRAMES [expr "$TOTAL_TRAJ_FRAMES - $prev_TOTAL_TRAJ_FRAMES"]
	    set prev_TOTAL_TRAJ_FRAMES $TOTAL_TRAJ_FRAMES
	    puts ""
	    puts "TRAJECTORY file ($c) loaded succesfully"
	    puts "  FILE:       $colvar_traj_file"
	    puts "  FRAMES:     $nFRAMES"
	    # check wether the number of frames in COLVAR file is synchronized with the number of frames in the trajectory XTC file
	    if {$COLVAR($c,traj_frames) != $nFRAMES} {
		puts "WARNING! the TRAJECTORY file is not synchronized with the corresponding COLVAR file"
		puts "         (number of loaded FRAMES should be $COLVAR($c,traj_frames)"
		puts "         it is very risky to continue with the analysis"
	    }
	}
	puts ""
	puts "TOTAL TRAJECTORY frames: $TOTAL_TRAJ_FRAMES"
	# check wether the total number of frames in all COLVAR file is synchronized with the total number of frames in all the trajectory XTC files
	if {$TOTAL_TRAJ_FRAMES != $TOTAL_COLVAR_STRUCTURES} {
	    puts "WARNING! total trajectory FRAMES ($TOTAL_TRAJ_FRAMES) differs from total colvar STRUCTURES ($TOTAL_COLVAR_STRUCTURES)"
	    puts "         are you sure TRAJECTORY files are synchronized with the corresponding COLVAR files?"
	    puts ""
	    puts "SOLVE THIS FIRST AND TRY AGAIN LATER!"
	}
	puts ""
    }

    #
    # compute microstates
    proc do_clusters {} {
	variable nTRAJECTORIES
	variable COLVAR
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable CV_INFO_GRID
	variable START_FOLDER
	variable CLUSTER_TYPE
	variable input_external
	variable cut_off
	variable sieving
	variable sievlog
	variable num_cluster
	variable efdim

	cd $START_FOLDER

	# get ative CVs from the user interface
	compute_active
	list_clusters 

	if {$nCV_ACTIVE==0} {
	    error "At least one CV should be marked as active ('Use' column)"
	}

	# set clusters dimensions
	set CLUST_SPACE ""
	for { set i 0} {$i < $nCV_ACTIVE } { incr i } {
	    set idum $CV_ACTIVE($i)
	    set_dimensions $idum $CV_INFO_GRID($idum,min) $CV_INFO_GRID($idum,max) $CV_INFO_GRID($idum,grid)
	    lappend CLUST_SPACE $idum
	}
	set_cluster_space $CLUST_SPACE

	# get active CVs according to the HILLS loaded
	array unset HILLS
	for { set h 0} {$h < [nHILLS] } { incr h } {
	    parse_hills $h
	}
	set_use_hills_default
	compute_active_hills

	if { $CLUSTER_TYPE == "grid" } {
	    set input_external GRID
	    set sievlog "0"
	    run_micro
	    write_frames "FRAMES"
	}
	if { $CLUSTER_TYPE == "read" } {
	    set sievlog "0"
	    read_data_from_cluster "MICRO_RUN"
	    write_frames "FRAMES"
	}
	if { $CLUSTER_TYPE == "external" } {
	    run_micro
	    write_frames "FRAMES"
	}
	if { $CLUSTER_TYPE == "kmed" } {
	    set input_external "KMEDOIDS $num_cluster"
	    run_micro
	    write_frames "FRAMES"
	}
	if { $CLUSTER_TYPE == "grmc" } {
	    set input_external "GROMACS $cut_off"
	    run_micro
	    write_frames "FRAMES"
	}
	if { $CLUSTER_TYPE == "dkmed" } {
	    set input_external "DENSITY  $cut_off $efdim"
	    run_micro
	    write_frames "FRAMES"
	}


	# write microstates in the user interface (MICROSTATES box-list)
	clusters_datalist_fill
    }

    #
    # read COLVAR files into memory. See comments to read_colvars. NOT
    # called for runtime variables, which are cleared.
    proc load_data {} {
	variable nTRAJECTORIES
	variable nRUNTIME
	variable COLVAR
	variable TRAJ_INFO
	variable TOTAL_COLVAR_FRAMES
	set TOTAL_COLVAR_FRAMES 0
	variable TOTAL_COLVAR_STRUCTURES
	set TOTAL_COLVAR_STRUCTURES 0
	variable TRAJ_SKIP
	variable CV_INFO_GRID
	variable START_FOLDER
	cd $START_FOLDER

	variable nVARIABLES
	array unset COLVAR

	variable in1
	variable out2

	puts "INFO: in load_data"
	
	# get the number of collective variables in COLVAR files
	set varlist [parse_colvars_metadata $TRAJ_INFO(0,colvar_file)]
	set nVARIABLES [llength $varlist]
	set nRUNTIME 0
	
	for {set i 1} {$i <= $nVARIABLES} {incr i} {
	    set CV_INFO_GRID($i,type) [lindex $varlist [expr {$i-1}]]
	    set COLVAR($i,min) 9999999
	    set COLVAR($i,max) -9999999
	}

	# read and load into memory each colvar file
	for { set c 0} {$c < $nTRAJECTORIES } { incr c } {
	    puts ""
	    puts "reading COLVAR file $TRAJ_INFO($c,colvar_file)" 
	    puts ""
	    read_active_colvars_header $c
	    read_colvars $c
	    set TOTAL_COLVAR_FRAMES [expr "$TOTAL_COLVAR_FRAMES + $COLVAR($c,nframes)"]
	    set TOTAL_COLVAR_STRUCTURES [expr "$TOTAL_COLVAR_STRUCTURES + $COLVAR($c,traj_frames)"]
	}
	puts ""
	puts "TOTAL COLVAR frames: $TOTAL_COLVAR_FRAMES"
	puts "TOTAL COLVAR structures: $TOTAL_COLVAR_STRUCTURES"
	puts ""

	# update user interface with information read from colvar files
	list_cvs $in1.cv.list
	list_cvs_vis $out2.cv.list
	#XXX
	variable out1
	list_hills_conv .metagui3.hlf.nb.analys.nb.wham.conv.list
	list_clusters
	# automatic build the default difusion matrix
	build_DD
    }

    #
    # build a default difusion matrix
    proc build_DD {} {
	variable nVARIABLES
	variable CV_INFO_GRID
	variable DD

	for { set i  1} {$i <= $nVARIABLES} { incr i } {
	    set dds [expr "( $CV_INFO_GRID($i,max) - $CV_INFO_GRID($i,min) ) / $CV_INFO_GRID($i,grid)"]
	    for { set j  1} {$j <= $nVARIABLES} { incr j } {
		set DD0 0
		if {$i==$j} {
		    set DD0 [expr "$dds*$dds"]
		}
		if { ![info exists DD($i,$j)] } {
		    set DD($i,$j)  $DD0
		}
	    }
	}
    }

    #
    # analyze the structure of COLVAR files. Sanity check. Return a
    # list of their names (possibly spaces), whose length can be used
    # as nVARIABLES. Does not change "global" variables.
    proc parse_colvars_metadata {colvarfile} {

	set ret {}
	set nV -1

	# get the number of collective variables in the first COLVAR
	# file (all colvar files are suposed to have the same
	# structure).
	set in [open $colvarfile r]
	while {[gets $in line] != -1} {
	    if { [lindex $line 0] == "#!" && [lindex $line 1] == "FIELDS" } {
		set nV 0
		set cvs [lreplace $line 0 1]; # drop #! and FIELDS
		foreach cv $cvs {
		    if { $cv == "time" || $cv == "vbias" } {
			continue
		    }
		    incr nV
		    # If name is "informative", remember it (?)
		    if { ! [regexp {^cv} $cv] && ! [regexp {^@} $cv] } {
			lappend ret $cv
			# set CV_INFO_GRID($nV,type) $cv
		    } else {
			lappend ret "unk"
		    }
		}
		break
	    }
	}
	close $in

	if {$nV == -1  } {
	    error "ERROR: likely missing FIELDS header in COLVAR file"
	} else {
	    puts "INFO: parse_colvar_metadata returning $ret"
	}

	return $ret
    }


    #
    # analyze the structure of HILLS files
    proc parse_hills {HILLS_index} {
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable START_FOLDER
	variable TRAJ_INFO
	cd $START_FOLDER

	set h $HILLS_index

	set in [open $TRAJ_INFO([hills2traj $h],hills_file) r]
	gets $in line

	# parsing header. get the number and id of active CVs on each HILLS file
	if { [lindex $line 0] == "#!" && [lindex $line 1] == "ACTIVE" } {
	    set dim  [lindex $line 2]
	    set TRAJ_INFO([hills2traj $h],hills_dim) $dim
	    set TRAJ_INFO([hills2traj $h],hills_cvs) ""
	    for {set i 0} {$i < $dim} {incr i} {
		set ia [lindex $line [expr "3 + $i"]]
		set ia [expr "0 + $ia "]; # tg ??
		lappend TRAJ_INFO([hills2traj $h],hills_cvs) $ia
	    }
	    # get the label of this HILLS file
	    set TRAJ_INFO([hills2traj $h],hills_label) [lindex $line end]
	} else {
	    puts "WARNING!  header line not found in HILLS file  $TRAJ_INFO([hills2traj $h],hills_file)"
	    puts "example of header line: #! ACTIVE 2  1 2  ABC"
	    puts "are you using a Plumed version older than 1.3 ???"
	    puts "please fix this problem and try again"
	}

	close $in

	puts ""
	puts "HILLS file ($HILLS_index) header read succesfully"
	puts "  FILE:       $TRAJ_INFO([hills2traj $h],hills_file)"
	puts "  LABEL:      $TRAJ_INFO([hills2traj $h],hills_label)"
	puts "  DIMENSIONS: $TRAJ_INFO([hills2traj $h],hills_dim)"
	puts "  CVs ID:     $TRAJ_INFO([hills2traj $h],hills_cvs)"
	puts ""
    }

    #
    # set which HILLS files to use (according to the CVs checked as active in the user interface)
    proc set_use_hills_default {} {
	variable TRAJ_INFO
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable START_FOLDER
	cd $START_FOLDER

	for { set h 0} {$h < [nHILLS] } { incr h } {
	    set goodall 1
	    foreach ia $TRAJ_INFO([hills2traj $h],hills_cvs) {
		set good 0
		for { set i 0} {$i < $nCV_ACTIVE } { incr i } {
		    if { $ia == $CV_ACTIVE($i) } {
			set good 1
		    }
		}
		set goodall [ expr " $goodall * $good" ]
	    }
	    set TRAJ_INFO([hills2traj $h],hills_use) $goodall
	    if {$TRAJ_INFO([hills2traj $h],hills_use)} {
		puts "HILLS file $TRAJ_INFO([hills2traj $h],hills_file) WILL BE USED FOR WHAM ANALYSIS"
	    }
	}
    }

    # read the ACTIVE lines, including bias label (currently unused)
    # and the list of active CVs to be shown in the GUI. This is only
    # done if "biased" is true for that trajectory. Fills
    # TRAJ_INFO(*,hills_cvs) and TRAJ_INFO(*,bias_label)
    proc read_active_colvars_header {c} {
	variable TRAJ_INFO
	puts "INFO: in read_active_colvars_header with arg $c"

	set biased $TRAJ_INFO($c,bias)
	
	set colvar_file $TRAJ_INFO($c,colvar_file)
	set found 0
	set in [open $colvar_file r]
	while {[gets $in line] != -1} {
	    if [ regexp {^#! ACTIVE} $line ] {
		set found 1
		break
	    }
	}
	close $in

	# ACTIVE mandatory for biased case
	if {$biased && ! $found } {
	    error       "ERROR! bias not defined for COLVAR file: $colvar_file
			 Header line not found.
			 Example of header line: #! ACTIVE 2  1 2  ABC
			 Plumed version >= 1.3 is required.
			 Please fix this problem and try again." 
	}

	# Fallbacks only for the unbiased case
	set bias_label "Unlabeled"
	set hills_cvs ""

	if {$found} {
	    # Now $line is as follows
	    # #! ACTIVE <N> <cv1> ... <cvN> <label>
	    # 0  1      2   3         3+N-1 end
	    set bias_label [lindex $line end]
	    set N [lindex $line 2]
	    set hills_cvs [lrange $line 3 [expr 2+$N]]
	}

	set TRAJ_INFO($c,bias_label) $bias_label
	set TRAJ_INFO($c,hills_cvs)  $hills_cvs
	
    }


    # Read and load into memory the contents of each COLVAR file.  The
    # argument $tr is a TRAJECTORY (!)  index, 0 < $tr < nTRAJECTORIES
    #
    # If the "runtime" arg is 1, the runtime colvars (there are
    # nRUNTIME of them) are appended to the existing ones
    # (1..nVARIABLES). The caller will later adjust nVARIABLES to
    # nVARIABLES+nRUNTIME . Even if runtime==1, the time column from
    # previously-loaded colvar files is used; this is due to the fact
    # that "driver" has no knowledge of real simulation timestamps.
    #
    # TRAJ_SKIP_BAFTI appears to be a way to stride frames, but it is
    # not bound to UI elements.
    #
    #  Reads:
    #    $TRAJ_INFO($tr,colvar_file)
    #    $TRAJ_INFO($tr,bias_label)
    #    possibly $TRAJ_INFO($tr,runtime_colvar)
    #  Sets:
    #    $COLVAR($tr,$f,time)
    #    $COLVAR($tr,$f,$i)    1 <= i <= nVARIABLES
    #    $COLVAR($i,min)        <--- note that the index is not tr!
    #    $COLVAR($i,max)        <--- note that the index is not tr!
    #    $COLVAR($tr,$f,bias_label)
    #    $COLVAR($tr,$f,cluster) -999
    #    $COLVAR($tr,nframes) $f
    #    $COLVAR($tr,traj_frames) $ff

    proc read_colvars {tr {runtime 0}} {
	variable COLVAR
	variable nVARIABLES
	variable nRUNTIME
	variable TRAJ_SKIP
	variable TRAJ_SKIP_BAFTI
	variable TRAJ_INFO
	variable T_CLUSTER
	variable START_FOLDER

	puts "INFO: in read_colvars with arguments $tr $runtime"
	
	cd $START_FOLDER

	set colvar_file $TRAJ_INFO($tr,colvar_file)
	set current_bias $TRAJ_INFO($tr,bias_label); # copied into COLVAR($tr,$f,bias)

	# Counters--
	set ii 0;		# actual lines in file
	set f 0;		# frames > T_CLUSTER
	set ff 0;		# frames > T_CLUSTER not skipped
	
	set in [open $colvar_file r]

	if {$runtime==0} {
	    # Read columns 1..nVARIABLES into CVs 1..nVARIABLES
	    set colvarOffset 0
	    set ncols $nVARIABLES
	} else {
	    # Read columns 1..nRUNTIME into CVs nVARIABLES+1 .. nVARIABLES+nRUNTIME
	    set colvarOffset $nVARIABLES
	    set ncols $nRUNTIME
	    
	    # If runtime, open and skip over comments
	    set in_r [open $TRAJ_INFO($tr,runtime_colvar) r]
	    while { [gets $in_r line_r] != -1} {
		if { ! [ regexp {^#} $line_r ] } {
		    break
		}
	    }
	}

	# Scan the rest of the COLVAR file
	while {[gets $in line] != -1} {
	    if [ regexp {^#} $line ] {
		# comment or header, do nothing
		continue
	    }

	    if {$runtime==0} {
		set line_data $line
	    } else {
		set line_data $line_r
		gets $in_r line_r
	    }

	    # Time comes from the "real" COLVAR file in any case
	    set this_time [lindex $line 0]

	    # read only COLVAR frames after T_CLUSTER
	    if {$this_time >= $T_CLUSTER} {
		# TONI I don't fully get the logic of first_traj_frame_to_be_loaded
		if {[info exists COLVAR($tr,first_traj_frame_to_be_loaded)]} {
		    # puts "DEBUG: first_traj_frame_to_be_loaded already set on trj. $tr"
		} else {
		    if {[expr "$ii%$TRAJ_SKIP_BAFTI"]==0} {
			set COLVAR($tr,first_traj_frame_to_be_loaded) [expr "$ii/$TRAJ_SKIP_BAFTI"]
		    }
		}
		
		if {[info exists COLVAR($tr,first_traj_frame_to_be_loaded)]} {
		    # These properties are being set again to their old value when runtime==1
		    set COLVAR($tr,$f,time) $this_time
		    set COLVAR($tr,$f,cluster) -999
			
		    for {set col 1} {$col <= $ncols } {incr col} {
			set cval [lindex $line_data $col]
			# i is the actual CV number, from 1
			set i [expr {$col+$colvarOffset}]
			set COLVAR($tr,$f,$i) $cval
			if {$cval < $COLVAR($i,min)} { set COLVAR($i,min) $cval }
			if {$cval > $COLVAR($i,max)} { set COLVAR($i,max) $cval }
			set COLVAR($tr,$f,bias) $current_bias; # TONI - I think this is redundant bc it's constant for each trajectory?
		    }
		    
		    incr f
		    if {[expr "($ii-$COLVAR($tr,first_traj_frame_to_be_loaded))%($TRAJ_SKIP*$TRAJ_SKIP_BAFTI)"] == 0} {
			incr ff
		    }
		}
	    }
	    incr ii
	}
    
	close $in
	if {$runtime==1} {
	    close $in_r
	}

	set COLVAR($tr,nframes) $f
	set COLVAR($tr,traj_frames) $ff
	puts ""
	puts "COLVAR file ($tr) read succesfully"
	puts "  FILE:        $colvar_file"
	puts "  VARIABLES:   $nVARIABLES"
	puts "  FRAMES:      $COLVAR($tr,nframes)"
	puts "  STRUCTURES:  $COLVAR($tr,traj_frames)"
	puts "  RUNTIME:     $runtime"
	if {$runtime==1} {
	    puts "    RUNTIME CVs: $nRUNTIME"
	    puts "    FILE:        $TRAJ_INFO($tr,runtime_colvar)"
	}
	puts ""
    }

    #
    # set space dimensions for each CV for microstates analysis
    proc set_dimensions {VARIABLE_number min_LIMIT max_LIMIT GRID_size} {
	variable DIMENSIONS

	set DIMENSIONS($VARIABLE_number,min) $min_LIMIT 
	set DIMENSIONS($VARIABLE_number,max) $max_LIMIT 
	set DIMENSIONS($VARIABLE_number,grid) $GRID_size 
	set DIMENSIONS($VARIABLE_number,spacer) [expr "(double($max_LIMIT)-double($min_LIMIT))/double($GRID_size)"]

	puts ""
	puts "GRID DIMENSIONS set for variable number $VARIABLE_number -> $DIMENSIONS($VARIABLE_number,min):$DIMENSIONS($VARIABLE_number,spacer):$DIMENSIONS($VARIABLE_number,max)"
	puts ""
    }

    #
    # set microstates space for each active CV
    proc set_cluster_space {VARIABLES_list} {
	variable CLUSTER_SPACE
	variable CLUSTER_SPACE_dim
	variable CLUSTER_SPACE_vol
	variable DIMENSIONS
	variable CV_INFO_GRID
	variable myVARS

	set CLUSTER_SPACE $VARIABLES_list
	set CLUSTER_SPACE_dim [llength $CLUSTER_SPACE]

	for {set i 0} {$i < $CLUSTER_SPACE_dim} {incr i} {
	    set myVARS($i,number) [lindex $CLUSTER_SPACE $i]
	    set myVARS($i,min)    $DIMENSIONS($myVARS($i,number),min)
	    set myVARS($i,max)    $DIMENSIONS($myVARS($i,number),max)
	    set myVARS($i,grid)   $DIMENSIONS($myVARS($i,number),grid)
	    set myVARS($i,spacer) $DIMENSIONS($myVARS($i,number),spacer)
	    # define period in myVARS to properly compute distances (ALEX)
	    set myVARS($i,period) 0
	    if {$CV_INFO_GRID($myVARS($i,number),periodic) == 1} {
		set myVARS($i,period) [expr "double($myVARS($i,max))-double($myVARS($i,min))"]
	    }
	}

	set j $CLUSTER_SPACE_dim
	set myshift 1
	set j [expr "$j-1"]
	set myVARS($j,shift) $myshift
	while {$j > 0} {
	    set myshift [expr "$myshift*$myVARS($j,grid)"]
	    set j [expr "$j-1"]
	    set myVARS($j,shift) $myshift
	}

	set l [expr "$CLUSTER_SPACE_dim - 1"]
	set CLUSTER_SPACE_vol [expr "$myVARS(0,shift)*$myVARS($l,grid)"]

	puts ""
	puts "CLUSTER SPACE set to variables numbers: $CLUSTER_SPACE"
	puts "CLUSTER SPACE dimensionality = $CLUSTER_SPACE_dim dimensions"
	puts "CLUSTER SPACE volume         = $CLUSTER_SPACE_vol clusters"
	puts ""
    }


    #
    # write input file microstates.in for clusters calculation
    # (MICROSTATES) with external program (cluster.x). Called from
    # run_micro
    proc write_cluster_input {{MICRO_file "microstates.in"}} {
	variable START_FOLDER
	variable CLUSTER_SPACE_dim
	variable nCV_ACTIVE
	variable CV_INFO_GRID
	variable CV_ACTIVE
	variable COLVAR
	variable myVARS
	variable nTRAJECTORIES
	variable START_FOLDER
	variable T_CLUSTER
	variable CLUSTER_TYPE
	variable input_external
	variable USE_ALL_WALKERS
	variable TRAJ_INFO
	
	variable cutting
	variable nconmin
	variable sievlog
	variable sieving

	cd $START_FOLDER
	cd MICRO_RUN


	# get active CVs according to the HILLS loaded
	array unset HILLS
	for { set h 0} {$h < [nHILLS] } { incr h } {
	    parse_hills $h
	}
	set_use_hills_default
	compute_active_hills
	set HILLS_in_use ""
	for { set h  0} {$h < [nHILLS] } { incr h } {
	    if { $TRAJ_INFO([hills2traj $h],hills_use) == 1 } {
		lappend HILLS_in_use $TRAJ_INFO([hills2traj $h],hills_label)
	    }
	}

	# generate clustering program input
	cd $START_FOLDER
	cd MICRO_RUN

	set out [open $MICRO_file w]

	puts $out "ARGUMENTS $input_external"
	puts $out "NCV       $CLUSTER_SPACE_dim"
	puts $out "CONNECTIVITY $cutting $nconmin "
	if {${sievlog}==1} {
	    puts $out "SIEVING $sieving"
	}
	set j 0
	for {set i 0} {$i < $CLUSTER_SPACE_dim} {incr i} {
	    incr j
	    puts $out  "CVGRID    $j $myVARS($i,min) $myVARS($i,max) $myVARS($i,grid)"
	}

	set PERCVS ""
	for { set i 0} {$i < $nCV_ACTIVE } { incr i } {
	    set idum $CV_ACTIVE($i)
	    if {$CV_INFO_GRID($idum,periodic)==0} {
		set PERCVS "$PERCVS .false."
	    }
	    if {$CV_INFO_GRID($idum,periodic)==1} {
		set PERCVS "$PERCVS .true."
	    }
	}
	puts $out "PERIODIC  $PERCVS"

	# write CV coordinates at times >= T_CLUSTER from all COLVAR files toguether
	set out2 [open "COORDINATES" w]
	set ff 0
	for {set c 0} {$c < $nTRAJECTORIES} {incr c} {
	    for {set f 0} {$f < $COLVAR($c,nframes)} {incr f} {
		incr ff
		# only assign frames biased by HILLS loaded and in use
		set good 1
        # 
        # Comment from Alex to discuss with Alessandro:
        # 
        # In bias-exchange trajectories, sometime ago, we arrived to the conclusion that we should only employ
        # the trajectories from the biased trajectories "in use" for generating the microstates. This makes,
        # of course, the non-biased trajectories impossible to analyze. Now I'm commenting that control until
        # we decide how to proceed. (set bb 1 per default) 
                set bb 1
	#	set bb 0
        # 
	#foreach h $HILLS_in_use {
	#    if { $COLVAR($c,$f,bias) == $h} { set bb 1}
	#}

		# forget it and use all walkers for the clustering
		if { $USE_ALL_WALKERS == 1 } { set bb 1 }
		
		# discard frames at times bellow T_CLUSTER
		if { $COLVAR($c,$f,time) < $T_CLUSTER } { set good 0 }
		# add label to frames outside grid boundaries
		for {set i 0} {$i < $CLUSTER_SPACE_dim} {incr i} {
		    if { $COLVAR($c,$f,$myVARS($i,number)) < $myVARS($i,min) } { set bb 0 }
		    if { $COLVAR($c,$f,$myVARS($i,number)) > $myVARS($i,max) } { set bb 0 }
		}

		if { $good == 1 } {
		    set CVs ""
		    for {set i 0} {$i < $CLUSTER_SPACE_dim} {incr i} {
			lappend CVs $COLVAR($c,$f,$myVARS($i,number))
		    }
		    puts $out2 "$f \t $c \t $CVs \t $COLVAR($c,$f,time) \t $bb"
		}
	    }
	}

	close $out2
	close $out
    }

    #
    # Link or, if impossible, copy, in a system-independent way. For
    # now just relies on TCL's file link. Link is deleted if
    # exists. Also handles the case when link is a directory, like "ln
    # -s" TG
    proc file_link_or_copy {src lnk} {
	if [file isdirectory $lnk] {
	    # {/tmp/src /dest} -> /dest/src
	    # {/tmp/src .} -> ./src
	    set lnk [file join $lnk [file tail $src]]
	}
	file delete $lnk;	# no problem if missing
	file link $lnk $src;	# note order
    }

    
    #
    # execute cluster.x command
    proc run_micro {{MICRO_file "microstates.in"} {MICRO_FOLDER "MICRO_RUN"}} {
	variable nTRAJECTORIES
	
	variable COLVAR
	
	variable START_FOLDER
	variable METAGUI_exe
	variable sievlog

	cd $START_FOLDER

	puts "cleaning $MICRO_FOLDER folder if exists"
	file delete -force $MICRO_FOLDER 
	#{*} [glob -nocomplain ./$MICRO_FOLDER/*]

	puts "moving to 'MICRO_RUN' folder"
        set ctrm [file exists $MICRO_FOLDER ]
        if { $ctrm == 0 } { file mkdir $MICRO_FOLDER }
	cd $MICRO_FOLDER

	puts "preparing clustering files"
	write_cluster_input $MICRO_file

	puts "running $METAGUI_exe  MICROSTATES $MICRO_file > microstates.out"
	cd $START_FOLDER
	cd $MICRO_FOLDER
	exec $METAGUI_exe MICROSTATES $MICRO_file > microstates.out
	puts "done"
	cd ..
	puts ""

	read_data_from_cluster
    }

    #
    # retrieve the results from a cluster.x calculation
    proc read_data_from_cluster {{MICRO_FOLDER "MICRO_RUN"}} {
	variable START_FOLDER
	cd $START_FOLDER

	puts "retrieving data from $MICRO_FOLDER"

	# TG: Link instead?
	file copy -force $MICRO_FOLDER/MICROSTATES .
	file copy -force $MICRO_FOLDER/connectivity .
	read_clusters_and_fes_and_basins "MICROSTATES"

	read_frames "$MICRO_FOLDER/FRAMES"

	puts "done"
	puts ""
    }

    #
    # write MICROSTATES file
    proc write_clusters {{CLUSTERS_file "CLUSTERS"}} {
	variable CLUSTERS
	variable NCL
	variable START_FOLDER
	cd $START_FOLDER
	
	set out [open $CLUSTERS_file w]

	for {set i 1} {$i <= $NCL} {incr i} {
	    puts -nonewline $out "$i \t $CLUSTERS($i,size) "
	    foreach cc $CLUSTERS($i,center) {
		puts -nonewline $out "\t $cc "
	    }
	    puts $out ""
	}
	close $out

	puts ""
	puts "$NCL CLUSTERS written to file $CLUSTERS_file"
	puts ""
    }

    #
    # write LABELS file
    proc write_labels {{LABELS_file "LABELS"}} {
	variable nTRAJECTORIES
	variable COLVAR
	variable T_CLUSTER
	variable START_FOLDER
	variable CLUSTER_SPACE_dim
	variable CV_ACTIVE_HILLS

	set owd [pwd]
	cd $START_FOLDER
	set out [open $LABELS_file w]

	set ff 0
	for {set c 0} {$c < $nTRAJECTORIES} {incr c} {
	    for {set f 0} {$f < $COLVAR($c,nframes)} {incr f} {
		if { $COLVAR($c,$f,time) >= $T_CLUSTER } {
		    if {$COLVAR($c,$f,cluster)>0} {
			set CVs ""
			foreach ia $CV_ACTIVE_HILLS {
			    lappend CVs $COLVAR($c,$f,$ia)
			}
			puts $out "$COLVAR($c,$f,time) \t $COLVAR($c,$f,cluster) \t $c \t $COLVAR($c,$f,bias) \t $CVs"
			incr ff
		    }
		}
	    }
	}

	close $out
	cd $owd

	puts ""
	puts "$ff frames written to file $LABELS_file"
	puts ""
    }

    #
    # wirte LABELS file including also frames before T_CLUSTER
    proc write_full_labels {{LABELS_file "LABELS_FULL"}} {
	variable nTRAJECTORIES
	variable COLVAR
	variable T_CLUSTER
	variable START_FOLDER
	cd $START_FOLDER

	set out [open $LABELS_file w]

	set ff 0
	for {set c 0} {$c < $nTRAJECTORIES} {incr c} {
	    for {set f 0} {$f < $COLVAR($c,nframes)} {incr f} {
		puts $out "$COLVAR($c,$f,time) \t $COLVAR($c,$f,cluster) \t $c \t $COLVAR($c,$f,bias)"
	    }
	}

	close $out

	puts ""
	puts "$ff frames written to file $LABELS_file"
	puts ""
    }

    #
    # write FRAMES file containing the microstate assignation at each loaded frame
    proc write_frames {{FRAMES_file "FRAMES"}} {
	variable FRAMES
	variable TRAJ_molid
	variable TOTAL_TRAJ_FRAMES
	variable START_FOLDER
	cd $START_FOLDER

	set out [open $FRAMES_file w]

	for {set f 0} {$f < $TOTAL_TRAJ_FRAMES} {incr f} {
	    puts $out "$f \t $FRAMES($f,cluster) \t $FRAMES($f,colvar_index) \t $FRAMES($f,colvar_frame) \t $FRAMES($f,colvar_time)"
	}

	close $out

	puts ""
	puts "$TOTAL_TRAJ_FRAMES frames written to $FRAMES_file"
	puts ""
    }

    #
    # write input for WHAM calculation according to user choices in the Graphical User Interface
    proc write_wham_input {{WHAM_file "wham.in"}} {
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable T_FILL
	variable DELTA
	variable KT
	
	
	
	variable nTRAJECTORIES
	variable COLVAR
	variable TRAJ_INFO
	variable T_CLUSTER
	variable DIMENSIONS
	variable CV_INFO_GRID
	variable myVARS
	variable START_FOLDER
	variable CV_ACTIVE_HILLS
	variable NCV_ACTIVE_HILLS
	variable GCORR
	variable TR_N_EXP  
	variable IGNORE_CONTROLS

	set out [open $WHAM_file w]

	puts $out "KT \t $KT"
	puts $out "T_FILL  \t $T_CLUSTER \t $T_FILL"
	puts $out "DELTA \t $DELTA"
	puts $out "GCORR \t $GCORR"
	puts $out "TR_N_EXP \t $TR_N_EXP"

	for { set h  0} {$h < [nHILLS] } { incr h } {
	    if { $TRAJ_INFO([hills2traj $h],hills_use) == 1 } {
		set tt ""
		foreach ia $TRAJ_INFO([hills2traj $h],hills_cvs) {
		    set jj 0
		    foreach ib $CV_ACTIVE_HILLS {
			incr jj
			if {  $ia==$ib } {
			    lappend tt   $jj 
			}
		    }
		}
		# TONI necessary bc. WHAM can't read slashes. I guess
		# this is the reason why symlinks or copies are done
		# in a WHAM_RUN subdir
		set bn [file tail $TRAJ_INFO([hills2traj $h],hills_file)]
		puts $out  "FILE_HILL \t $bn \t $TRAJ_INFO([hills2traj $h],hills_dim) \t $tt"
	    }
	}

	puts $out "NCV_ACTIVE \t $NCV_ACTIVE_HILLS"


	set j 0
	foreach cv $CV_ACTIVE_HILLS {
	    incr j
	    puts $out "RANGE_CV $j $CV_INFO_GRID($cv,min) $CV_INFO_GRID($cv,max)"
	    if {$CV_INFO_GRID($cv,periodic)==1} { 
		puts $out "PERIODIC $j" 
	    }
	}

	puts $out "NGRID \t 50"
	puts $out "MIN_OVERLAP 1"
	if {$IGNORE_CONTROLS==1} {
	    puts $out "IGNORE_CONTROLS"
	}
	
	close $out
    }

    #
    # execute wham.x command
    proc run_wham {{WHAM_FOLDER_rel "WHAM_RUN"}} {
	variable nTRAJECTORIES
	variable TRAJ_INFO
	variable COLVAR
	variable START_FOLDER
	variable METAGUI_exe
	variable nCV_ACTIVE
	variable winf

	if {$nCV_ACTIVE==0} {
	    error "At least one CV should be marked as active ('Use' column)"
	}

	set WHAM_FOLDER [file join $START_FOLDER $WHAM_FOLDER_rel]

	cd $START_FOLDER
	puts "cleaning $WHAM_FOLDER folder if exists"
	file delete -force $WHAM_FOLDER; # does not fail
	file mkdir $WHAM_FOLDER;	 # may fail (e.g. permissions)
	cd $WHAM_FOLDER

	puts "preparing wheighted-histogram analysis files"
	write_labels
	write_wham_input

	# "file join" needed to handle absolute pathnames
	for { set h 0} {$h < [nHILLS] } { incr h } {
	    file_link_or_copy [file join $START_FOLDER $TRAJ_INFO([hills2traj $h],hills_file)] .
	}

	# TONI are COLVARs needed for WHAM? disabling them
	#for { set c 0} {$c < $nTRAJECTORIES } { incr c } {
        #	    file_link_or_copy [file join $START_FOLDER $TRAJ_INFO($c,colvar_file)] .
	#}

	file copy -force ../MICROSTATES ./CLUSTERS; # TG link?
	file copy -force ../LABELS .;		    # TG link?
	puts "running $METAGUI_exe WHAM wham.in > wham.out"
	exec $METAGUI_exe WHAM wham.in > wham.out 
	puts "done"

	puts ""
	read_data_from_wham
	clusters_datalist_fill
    }

    #
    # retrieve the results from a wham.x calculation
    proc read_data_from_wham {{WHAM_FOLDER "WHAM_RUN"}} {
	variable START_FOLDER
	cd $START_FOLDER

	puts "retrieving data from $WHAM_FOLDER"
	read_fes "$WHAM_FOLDER/FES"
	write_clusters_and_fes_and_basins "MICROSTATES"
	put_fes_on_user
	puts ""
    }

    #
    # update MICROSTATES file with the free energy of each microstate
    proc write_clusters_and_fes {{CLUSTERS_and_FES_file "CLUSTERS.FES"}} {
	variable CLUSTERS
	variable NCL
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable DIMENSIONS
	variable START_FOLDER
	cd $START_FOLDER

	set out [open $CLUSTERS_and_FES_file w]

	for {set i 1} {$i <= $NCL} {incr i} {
	    puts -nonewline $out "$i \t $CLUSTERS($i,size) "
	    foreach cc $CLUSTERS($i,center) {
		puts -nonewline $out "\t $cc "
	    }
	    puts -nonewline $out "\t $CLUSTERS($i,fes) "
	    puts $out ""
	}
	close $out

	puts ""
	puts "$NCL CLUSTERS written to file $CLUSTERS_and_FES_file"
	puts ""
    }

    #
    # write input files for SAF cluster
    proc write_saf_input {{CFILE "SAF.in"}} {
	variable CLUSTERS
	variable NCL
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable DIMENSIONS
	variable START_FOLDER
	variable KT
	variable myVARS

	set out [open $CFILE w]
	puts $out "MICROSTATES"
	puts $out "$KT"
	set COORDINATES_dim [llength $CLUSTERS(1,center)]
	puts $out "$COORDINATES_dim"
	for {set i 0} {$i < $COORDINATES_dim} {incr i 1} {
	    set myspacer $myVARS($i,spacer)
	    set myperiod $myVARS($i,period)
	    set j [expr "$i + 1"]
	    puts $out "$j $myperiod $myspacer"
	}

	close $out
    }


    proc generate_input_smart {{SAF_FOLDER "SAF_RUN"}} {
	# variables for smart clustering
	variable RHO_CUT
	variable DELTA_CUT
	variable START_FOLDER
	cd $START_FOLDER
	set in [open "$START_FOLDER/$SAF_FOLDER/tmp" w ]
	puts $in "$RHO_CUT"
	puts $in "$DELTA_CUT"
	puts $in " "
	close $in
	update
	# TG link?
	file copy -force $START_FOLDER/$SAF_FOLDER/tmp $START_FOLDER/$SAF_FOLDER/input.smart
	# infinite loop waiting for $SAF_exe to finish
	set contr 0
	while {$contr<1} {
	    set contr [file exists $SAF_FOLDER/smart.out]
	    after 1000
	}
	file delete -force $SAF_FOLDER/input_SAF
	build_basins_SAF
    }

    #
    # read output data from  smart clusters and assign basins
    proc build_basins_SAF {{FES_RANGE 100} {SAF_FOLDER "SAF_RUN"}} {
	variable NCL
	variable CLUSTERS
	variable KT
	variable NB
	variable BASINS
	variable nCV_ACTIVE
	variable START_FOLDER
	variable var_list

	cd $START_FOLDER

	array unset BASINS

	# read number of required basins subdivision
	set in [open "$SAF_FOLDER/smart.out" r]
	gets $in line
	set NB [lindex $line 0]

	# initialize variables
	for {set i 1} {$i <= $NCL} {incr i} {
	    set CLUSTERS($i,basin) -1000
	}
	for {set b 1} {$b <= $NB} {incr b} {
	    set BASINS($b,free_min) 9999999
	    set BASINS($b,avrg_fes) 0.000
	    set BASINS($b,clusters_list) ""
	    set BASINS($b,size) 0
	    set BASINS($b,cluster_center) 0
	}


	while {[gets $in line] != -1} {
	    set cl [lindex $line 0]
	    set b  [lindex $line 1]    
	    if {$b!=0} {
		set bas $b
		# attractor of each microstate is indexed as a negative number
		if {$b<0} {
		    set bas [expr "-1*$b"]	
		    set BASINS($bas,cluster_center) $cl
		    set BASINS($bas,free_min) $CLUSTERS($cl,fes)
		}
		# each microstate is assigned to a kinetic basin
		set CLUSTERS($cl,basin) $bas
	    }
	}
	close $in

	set total_ass 0
	set BASINS(-999,size) 0
	set BASINS(-999,clusters_list) ""

	# get the size (in number of microstates) of each kinetic basin and compute average free energies
	for {set i 1} {$i <= $NCL} {incr i} {
	    set free $CLUSTERS($i,fes)
	    set b $CLUSTERS($i,basin)
	    if {$b>0} {
		set free_min $BASINS($b,free_min)
		set free_diff [expr "$free - $free_min"]
		set free_max [expr "$FES_RANGE"]
		if { $free_diff < $free_max } {
		    incr BASINS($b,size)
		    set BASINS($b,avrg_fes) [expr "$BASINS($b,avrg_fes) + exp(-$free/$KT)"]
		    lappend BASINS($b,clusters_list) $i
		    incr total_ass
		} else {
		    set CLUSTERS($i,basin) -999
		    incr BASINS(-999,size)
		    lappend BASINS(-999,clusters_list) $i
		}
	    }
	}
	#
	puts "$NB BASINS found"
	puts "$total_ass CLUSTERS (out of $NCL) assigned to the BASINS"
	#
	for {set b 1} {$b <= $NB} {incr b} {
	    if {$BASINS($b,size)>0} {
		set BASINS($b,avrg_fes) [expr "-$KT*log($BASINS($b,avrg_fes)/$BASINS($b,size))"]
	    } else {
		puts "WARNING: basin number $b is empty"
	    }
	}
	# write basins assignation to BASINS file
	write_basins "BASINS"
	# update MICROSTATES file with basins assignation
	write_clusters_and_fes_and_basins "MICROSTATES"
	clusters_datalist_fill
	basins_datalist_fill
    }

    # RUNING SAF
    proc run_SAF {{SAF_FOLDER "SAF_RUN"}} {
	variable METAGUI_exe
	variable START_FOLDER
	cd $START_FOLDER
	puts "cleaning $SAF_FOLDER folder if exists"
	file delete -force $SAF_FOLDER 
	#file delete -force $START_FOLDER/$SAF_FOLDER/input.smart
	#file delete -force $START_FOLDER/$SAF_FOLDER/MICROSTATES

	puts "moving to '$SAF_FOLDER' folder"

        set ctrm [file exists $SAF_FOLDER]
	if { $ctrm == 0 } { file mkdir $SAF_FOLDER }
	file copy -force $START_FOLDER/MICROSTATES $START_FOLDER/$SAF_FOLDER/MICROSTATES
	cd $SAF_FOLDER
	write_saf_input
	exec $METAGUI_exe CLUSTERS SAF.in >SAF.out &
	set contr 0
	# infinite loop waiting dec.dat to be generated- Toni: why?
	while {$contr<1} {
	    set contr [file exists $START_FOLDER/$SAF_FOLDER/dec.dat]
	    after 1000
	}
	show_DEC
    }


    #
    # write input and execute the kinetic_basins.x command. this will perform a default subdivision in 1 to 15 kinetic basins
    proc run_basins {{CLUSTERS_and_FES_file "CLUSTERS.FES"} {MAX_REQUESTED_EIGENVECTORS 15} {BASINS_FOLDER "BASINS_RUN"}} {
	variable CLUSTER_SPACE_dim
	variable myVARS
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable KT
	variable START_FOLDER
	variable DD
	variable REFES
	variable METAGUI_exe
	variable CV_INFO_GRID   

	cd $START_FOLDER

	set nDIM $nCV_ACTIVE

	# preparing input file for kinetic_basins.x analysis according to user choices in the GUI

	set ds_string ""
	for {set c 0} {$c < $nDIM} {incr c} {
	    set ds_string "$ds_string $myVARS($c,spacer)"
	}

	set grid_string ""
	for {set c 0} {$c < $nDIM} {incr c} {
	    set grid_string "$grid_string $myVARS($c,grid)"
	}

	set j 0
	set PERCVS ""
	for {set c 0} {$c < $nDIM} {incr c} {
	    set idum $CV_ACTIVE($c)
	    if {$CV_INFO_GRID($idum,periodic)==0} {
		set PERCVS "$PERCVS .false."
	    }
	    if {$CV_INFO_GRID($idum,periodic)==1} {
		set PERCVS "$PERCVS .true."
	    }
	}

	set FESREF ""

	if {$REFES==1} {
	    set FESREF "$FESREF .true."
	}
	
	if {$REFES==0} {
	    set FESREF "$FESREF .false."
	}    

	puts "cleaning $BASINS_FOLDER folder if exists"
	file delete -force $BASINS_FOLDER
	puts "moving to '$BASINS_FOLDER' folder"
        set ctrm [file exists $BASINS_FOLDER ]
        if { $ctrm == 0 } { file mkdir $BASINS_FOLDER }
	# write kinetic_basins.x input file
	cd $BASINS_FOLDER

	puts "preparing kinetic matrix analysis files"

	set in [open "kin.in" w]
	puts $in "$CLUSTERS_and_FES_file"
	puts $in "$nDIM $MAX_REQUESTED_EIGENVECTORS $KT"
	puts $in "$ds_string"
	puts $in "$grid_string"
	puts $in "$PERCVS"

	for {set c 0} {$c < $nDIM} {incr c} {
	    set DDall ""
	    for {set j 0} {$j < $nDIM} {incr j} {
		set DDall "$DDall $DD($myVARS($c,number),$myVARS($j,number))"
	    }
	    puts $in "$DDall"
	}
	puts $in "$FESREF"    
	close $in
	file copy -force ../MICROSTATES ./CLUSTERS.FES
	file copy -force ../connectivity .
	# execute kinetic_basins.x
	puts "running $METAGUI_exe KINETIC kin.in > kin.out"
	exec $METAGUI_exe KINETIC kin.in > kin.out
	puts "done"
	cd ..
    }


    #
    # read output data from  kinetic_basins.x analysis and assing microstates to each kinetic basins they belong to
    proc build_basins {{NUMBER_of_VECTORS 2} {FES_RANGE 100} {EIGENVECTORS_folder "BASINS_RUN"}} {
	variable NCL
	variable CLUSTERS
	variable KT
	variable NB
	variable BASINS
	variable nCV_ACTIVE
	variable START_FOLDER
	variable var_list

	cd $START_FOLDER

	array unset BASINS

	# read number of required basins subdivision
	set nVECT $NUMBER_of_VECTORS
	set NB [expr "$NUMBER_of_VECTORS + 1"]

	# initialize variables
	for {set i 1} {$i <= $NCL} {incr i} {
	    set CLUSTERS($i,basin) -1000
	}
	for {set b 1} {$b <= $NB} {incr b} {
	    set BASINS($b,free_min) 9999999
	    set BASINS($b,avrg_fes) 0.000
	    set BASINS($b,clusters_list) ""
	    set BASINS($b,size) 0
	    set BASINS($b,cluster_center) 0
	}

	puts ""
	puts "reading BASINS made of $NUMBER_of_VECTORS EIGENVECTORS from $EIGENVECTORS_folder/CLUSTERS.BASINS.ALL_EIGENVECTORS"
	puts ""

	set in [open "$EIGENVECTORS_folder/CLUSTERS.BASINS.ALL_EIGENVECTORS" r]

	# read basins assignation from kinetic_basins.x output file
	while {[gets $in line] != -1} {
	    set cl [lindex $line 0]
	    set b  [lindex $line $NUMBER_of_VECTORS]    
	    if {$b!=0} {
		set bas $b
		# attractor of each microstate is indexed as a negative number
		if {$b<0} {
		    set bas [expr "-1*$b"]	
		    set BASINS($bas,cluster_center) $cl
		    set BASINS($bas,free_min) $CLUSTERS($cl,fes)
		}
		# each microstate is assigned to a kinetic basin
		set CLUSTERS($cl,basin) $bas
	    }
	}
	close $in

	set total_ass 0
	set BASINS(-999,size) 0
	set BASINS(-999,clusters_list) ""

	# get the size (in number of microstates) of each kinetic basin and compute average free energies
	for {set i 1} {$i <= $NCL} {incr i} {
	    set free $CLUSTERS($i,fes)
	    set b $CLUSTERS($i,basin)
	    if {$b>0} {
		set free_min $BASINS($b,free_min)
		set free_diff [expr "$free - $free_min"]
		set free_max [expr "$FES_RANGE"]
		if { $free_diff < $free_max } {
		    incr BASINS($b,size)
		    set BASINS($b,avrg_fes) [expr "$BASINS($b,avrg_fes) + exp(-$free/$KT)"]
		    lappend BASINS($b,clusters_list) $i
		    incr total_ass
		} else {
		    set CLUSTERS($i,basin) -999
		    incr BASINS(-999,size)
		    lappend BASINS(-999,clusters_list) $i
		}
	    }
	}

	puts "$NB BASINS found"
	puts "$total_ass CLUSTERS (out of $NCL) assigned to the BASINS"

	for {set b 1} {$b <= $NB} {incr b} {
	    if {$BASINS($b,size)>0} {
		set BASINS($b,avrg_fes) [expr "-$KT*log($BASINS($b,avrg_fes)/$BASINS($b,size))"]
	    } else {
		puts "WARNING: basin number $b is empty"
	    }
	}
	
	# write basins assignation to BASINS file
	write_basins "BASINS"
	# update MICROSTATES file with basins assignation
	write_clusters_and_fes_and_basins "MICROSTATES"

	clusters_datalist_fill
	basins_datalist_fill
    }

    #
    # write basins assignation to BASINS file
    proc write_basins {{BASINS_file "BASINS"}} {
	variable CLUSTERS
	variable NCL
	variable KT
	variable NB
	variable BASINS
	variable START_FOLDER
	cd $START_FOLDER

	set out [open $BASINS_file w]

	for {set b 1} {$b <= $NB} {incr b} {
	    if {$BASINS($b,size)>0} {
		puts -nonewline $out "$b \t $BASINS($b,size) \t $BASINS($b,avrg_fes) "

		set cl $BASINS($b,cluster_center)
		puts -nonewline $out "\t $cl "
		foreach cc $CLUSTERS($cl,center) {
		    puts -nonewline $out "\t $cc "
		}
		puts -nonewline $out "\t $CLUSTERS($cl,fes) "
		puts $out ""
	    }
	}
	close $out

	puts ""
	puts "$NB BASINS written to file $BASINS_file"
	puts ""
    }

    #
    # update MICROSTATES file with basins assignation
    proc write_clusters_and_fes_and_basins {{CLUSTERS_and_FES_and_BASINS_file "CLUSTERS.FES.BASINS"}} {
	variable CLUSTERS
	variable NCL
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable DIMENSIONS
	variable START_FOLDER
	cd $START_FOLDER

	set out [open $CLUSTERS_and_FES_and_BASINS_file w]

	for {set i 1} {$i <= $NCL} {incr i} {
	    puts -nonewline $out "$i \t $CLUSTERS($i,size) "
	    foreach cc $CLUSTERS($i,center) {
		puts -nonewline $out "\t [format %.8f $cc] "
	    }
	    puts -nonewline $out "\t $CLUSTERS($i,fes) "
	    puts -nonewline $out "\t $CLUSTERS($i,basin) "
	    puts $out ""
	}
	close $out

	puts ""
	puts "$NCL CLUSTERS, FES and BASINS written to file $CLUSTERS_and_FES_and_BASINS_file"
	puts ""
    }


    #
    # write and load the microstates representations (spheres) into VMD display
    proc show_basins {{SHOW_ACTIVE_VARS "0 1 2"} {FES_RANGE "100"} {TASK "BASINS"}} {
	variable myVARS
	variable START_FOLDER
	cd $START_FOLDER
	if {$SHOW_ACTIVE_VARS=="error" || $SHOW_ACTIVE_VARS=="F F F"} {
	    error "The coordinates for plotting are not properly set"
	    # return "the coordinates for plotting are not properly set"
	}

	write_PDB_clusters_and_fes_and_basins "MICROSTATES.pdb" "$SHOW_ACTIVE_VARS" $TASK $FES_RANGE
	load_fes "MICROSTATES.pdb" $TASK $FES_RANGE
    }

    #
    # reads basins assignation from BASINS file. 
    # TODO: not currently assigned to any button in the GUI.
    # TODO: check consistency
    proc read_basins {{BASINS_file "BASINS"}} {
	variable BASINS
	variable NB
	variable START_FOLDER
	cd $START_FOLDER

	if {[info exists NB]} { set old_NB $NB } else { set old_NB 0 }

	puts ""
	puts "resetting old BASINS size, centers, fes and basins"
	for {set i 1} {$i <= $old_NCL} {incr i} {
	    if {[info exists CLUSTERS($i,size)]} {unset CLUSTERS($i,size)}
	    if {[info exists CLUSTERS($i,center)]} {unset CLUSTERS($i,center)}
	    if {[info exists CLUSTERS($i,fes)]} {unset CLUSTERS($i,fes)}
	    if {[info exists CLUSTERS($i,basin)]} {unset CLUSTERS($i,basin)}
	}
	puts ""


    }

    #
    # reads microstates data from MICROSTATES file. 
    proc read_clusters_and_fes_and_basins {{CLUSTERS_and_FES_and_BASINS_file "MICROSTATES"}} {
	variable CLUSTERS
	variable NCL
	variable nCV_ACTIVE
	variable FES_MIN
	variable FES_MAX
	variable START_FOLDER
	cd $START_FOLDER

	set nDIM $nCV_ACTIVE

	if {[info exists CLUSTERS]} { unset CLUSTERS }
	#    if {[info exists NCL]} { set old_NCL $NCL } else { set old_NCL 0 }
	#
	#    puts ""
	#    puts "resetting old CLUSTERS size, centers and fes"
	#    for {set i 1} {$i <= $old_NCL} {incr i} {
	#	if {[info exists CLUSTERS($i,size)]} {unset CLUSTERS($i,size)}
	#	if {[info exists CLUSTERS($i,center)]} {unset CLUSTERS($i,center)}
	#	if {[info exists CLUSTERS($i,fes)]} {unset CLUSTERS($i,fes)}
	#    }
	#    puts ""
	set CLUSTERS(-999,size) 0
	set CLUSTERS(-999,fes) -999
	set CLUSTERS(-999,bin) -999
	set CLUSTERS(-999,basin) -999
	set CLUSTERS(-999,frames) "" 

	set FES_MIN 999999999
	set FES_MAX -999999999

	set new_NCL 0
	set in [open $CLUSTERS_and_FES_and_BASINS_file r]
	while {[gets $in line] != -1} {
	    incr new_NCL
	    set clid [lindex $line 0]
	    set clsize [lindex $line 1]
	    set CLUSTERS($clid,size) $clsize
	    set CLUSTERS($clid,center) ""
	    for {set d 0} {$d < $nDIM} { incr d} { 
		set col [expr "2 + $d"] 
		lappend CLUSTERS($clid,center) [lindex $line $col]
	    }
	    set col [expr "2 + $nDIM"] 
	    set tmp [lindex $line $col]
	    set fes [expr "$tmp + 0.0"]
	    if {$fes < $FES_MIN} {set FES_MIN $fes}
	    if {$fes > $FES_MAX} {set FES_MAX $fes}
	    set CLUSTERS($clid,fes) $fes
	    set col [expr "2 + $nDIM + 1"] 
	    set basin [lindex $line $col]
	    set CLUSTERS($clid,basin) $basin
	    set CLUSTERS($clid,frames) " " 
	}

	puts ""
	puts "$new_NCL CLUSTERS, FES and BASINS read from $CLUSTERS_and_FES_and_BASINS_file"
	puts ""
	if { $new_NCL != $clid } { puts "WARNING! something went wrong. There are missing clusters ids in $CLUSTERS_and_FES_and_BASINS_file" }
	set NCL $new_NCL

    }

    #
    # reads frames assignation to microstates from FRAMES file.
    # requires read_clusters_and_fes_and_basins first.
    proc read_frames {{FRAMES_file "FRAMES"}} {
	variable CLUSTERS
	variable FRAMES
	variable TRAJ_molid
	variable TOTAL_TRAJ_FRAMES
	variable NCL
	variable START_FOLDER
	variable COLVAR
	variable nTRAJECTORIES
	variable TRAJ_SKIP
	variable TRAJ_SKIP_BAFTI

	cd $START_FOLDER

	puts ""
	puts "resetting frames assignment to CLUSTERS"
	puts ""
	if {[info exists FRAMES]} { unset FRAMES }

	# traj_ff is the vmd trajectory frame index (all TRAJECTORY files)
	set traj_ff 0
	for {set c 0} {$c < $nTRAJECTORIES} {incr c} {
	    for {set f 0} {$f < $COLVAR($c,nframes)} {incr f} {
		set COLVAR($c,$f,cluster) -999
		# this is special when performing the analysis on a subset of colvar frames (i.e. more COLVAR frames loaded than TRAJECTORY frames)
		set do_frame [expr "$f%($TRAJ_SKIP*$TRAJ_SKIP_BAFTI)"]
		if {$do_frame == 0} {
		    set FRAMES($traj_ff,cluster) -999
		    set FRAMES($traj_ff,colvar_index) $c
		    set FRAMES($traj_ff,colvar_frame) $f
		    set FRAMES($traj_ff,colvar_time) $COLVAR($c,$f,time)
		    incr traj_ff
		}
	    }
	}

	puts ""
	puts "reading frames from $FRAMES_file"
	puts ""
	set nFRAMES 0
	set traj_ff 0
	set in [open $FRAMES_file r]
	while {[gets $in line] != -1} {
	    incr nFRAMES
	    set ff [lindex $line 0]
	    set cl [lindex $line 1]
	    set colvar_index [lindex $line 2]
	    set colvar_frame [lindex $line 3]
	    set colvar_time [lindex $line 4]
	    set COLVAR($colvar_index,$colvar_frame,cluster) $cl
	    # this is special when performing the analysis on a subset of colvar frames (i.e. more COLVAR frames loaded than TRAJECTORY frames)
	    set do_frame [expr "$colvar_frame%($TRAJ_SKIP*$TRAJ_SKIP_BAFTI)"]
	    if {$do_frame == 0} {
		lappend CLUSTERS($cl,frames) $traj_ff
		set FRAMES($traj_ff,cluster) $cl
		set FRAMES($traj_ff,colvar_index) $colvar_index
		set FRAMES($traj_ff,colvar_frame) $colvar_frame
		set FRAMES($traj_ff,colvar_time) $colvar_time
		incr traj_ff
	    }
	}

	puts ""
	puts "$nFRAMES frames read from $FRAMES_file"
	if {$nFRAMES != $TOTAL_TRAJ_FRAMES} {
	    puts "WARNING! number of frames read ($traj_ff) differs from total trajectory frames ($TOTAL_TRAJ_FRAMES)"
	    puts "         this is not an error if you are not using all frames for the clustering."
	    puts "         FRAMES file will be written containing all loaded TRAJECTORY frames."
	}
	puts ""

    }

    #
    # reads MICROSTATES, BASINS and LABELS files. 
    # TODO: not currently assigned to any button in the GUI.
    # TODO: check consistency
    proc read_all_data {{DATA_FILE "CLUSTERS.FES"}} {
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable CV_INFO_GRID
	variable START_FOLDER
	cd $START_FOLDER

	set nDIM $nCV_ACTIVE

	for { set i 0} {$i < $nCV_ACTIVE } { incr i } {
	    set idum $CV_ACTIVE($i)
	    set_dimensions $idum $CV_INFO_GRID($idum,min) $CV_INFO_GRID($idum,max) $CV_INFO_GRID($idum,grid)
	    lappend CLUST_SPACE $idum
	}
	set_cluster_space $CLUST_SPACE

	read_clusters_fes $DATA_FILE
	read_frames
	put_fes_on_user
	puts ""
	puts "ALL DATA loaded succesfully from $DATA_FILE"
	puts "you should now re-build the basins (provided you computed the eigenvectors previously)"
	puts ""
    }

    #
    # write PDB file containing a representation of the microstates and basins subdivision of the system. 
    proc write_PDB_clusters_and_fes_and_basins {{CLUSTERS_and_FES_and_BASINS_PDBfile "CLUSTERS.FES.BASINS.pdb"} {WRITE_ACTIVE_VARS_list "0 1 2"} {TASK "BASINS"} {FES_RANGE "100"}} {
	variable CLUSTERS
	variable BASINS
	variable NCL
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable myVARS
	variable FES_MIN
	variable FES_MAX
	variable START_FOLDER
	cd $START_FOLDER

	set Xvar [lindex $WRITE_ACTIVE_VARS_list 0]
	set Yvar [lindex $WRITE_ACTIVE_VARS_list 1]
	set Zvar [lindex $WRITE_ACTIVE_VARS_list 2]

	set Fmin $FES_MIN
	
	if {$FES_RANGE=="100"} {
	    set Fmax $FES_MAX
	} else {
	    set Fmax $FES_RANGE
	}
	#BBBBB
	if {$TASK=="MICROSTATES"} {
	    set Fmax $FES_MAX
	} 

	set Xmin $myVARS($Xvar,min)
	set Xmax $myVARS($Xvar,max)
	if {$Yvar!="F"} {
	    set Ymin $myVARS($Yvar,min)
	    set Ymax $myVARS($Yvar,max)
	}
	if {$Zvar!="F" && $Zvar!="zero" } {
	    set Zmin $myVARS($Zvar,min)
	    set Zmax $myVARS($Zvar,max)
	}

	set out [open $CLUSTERS_and_FES_and_BASINS_PDBfile w]

	puts $out "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1"
	puts $out "REMARK     varX=$myVARS($Xvar,number)     minX=$Xmin   maxX=$Xmax"
	if {$Yvar=="F"} {
	    puts $out "REMARK     varY=F     minY=$Fmin   maxY=$Fmax"
	} else {
	    puts $out "REMARK     varY=$myVARS($Yvar,number)     minY=$Ymin   maxY=$Ymax"
	}
	if {$Zvar=="F"} {
	    puts $out "REMARK     varZ=F     minZ=$Fmin   maxZ=$Fmax"
	} elseif {$Zvar=="zero"} {
	    puts $out "REMARK     varZ=0     minZ=NA      maxZ=NA"
	} else {
	    puts $out "REMARK     varZ=$myVARS($Zvar,number)     minZ=$Zmin   maxZ=$Zmax"
	}
	puts $out "REMARK     F          minF=$Fmin   maxF=$Fmax"

	# TG: Note that Occ and Beta are misaligned. Does not seem to
	# cause problems.
	set PDB_MASK "%4s  %5i  %2s  %3s %5i    %8.3f%8.3f%8.3f  %3.2f  %3.2f         %2s"
	set CON_MASK "%6s%5i%5i%5i%5i%5i"
	#BBBBB
	# Number of atoms (microstates) multiplied by 5 in such a way that we can represent high connectivities

	set ifinta "0"
	for {set j 1} {$j <= 5 } {incr j} {
	    for {set i 1} {$i <= $NCL} {incr i} {
		#BBBBB
		set ifinta [expr "$ifinta+1"]
		if {$TASK!="MICROSTATES"} {
		    set F $CLUSTERS($i,fes)
		    set myBeta [expr "($F-$Fmin)/($Fmax-$Fmin)"]
		    if {$F<1000.} {
			set myResName "MIC"
		    } else {
			set myResName "BAD"
			set myBeta 1.
		    }
		} else {
		    set myResName "MIC"
		    set myBeta 1.
		}

		set X [lindex $CLUSTERS($i,center) $Xvar]
		set myX [expr "($X-$Xmin)/($Xmax-$Xmin)*100"]

		if {$Yvar=="F"} {
		    set myY [expr "$myBeta*100"]
		} else {
		    set Y [lindex $CLUSTERS($i,center) $Yvar]
		    set myY [expr "($Y-$Ymin)/($Ymax-$Ymin)*100"]
		}

		if {$Zvar=="F"} {
		    set myZ [expr "$myBeta*100"]
		} elseif {$Zvar=="zero"} {
		    set myZ 0.0
		} else {
		    set Z [lindex $CLUSTERS($i,center) $Zvar]
		    set myZ [expr "($Z-$Zmin)/($Zmax-$Zmin)*100"]
		}

		if {$TASK=="BASINS"} {
		    set B $CLUSTERS($i,basin)
		} else {
		    set B 1
		}
		set myOccup $CLUSTERS($i,fes)
		set myResid $B
		set myName "O"
		if {$TASK=="BASINS" && $B>0} {
		    if {$BASINS($B,cluster_center)==$i} {
			set myName "C"
		    }
		}
		if {$nCV_ACTIVE>3} {
		    set scaler 0
		    set my_rand [expr "5*(-1 + 2*rand())"]
		    set myX [expr "$my_rand * $scaler + $myX"]
		    if {$Yvar!="F"} {
			set my_rand [expr "5*(-1 + 2*rand())"]
			set myY [expr "$my_rand * $scaler + $myY"]
		    }
		    if {$Zvar!="F"} {
			set my_rand [expr "5*(-1 + 2*rand())"]
			set myZ [expr "$my_rand * $scaler + $myZ"]
		    }
		}

		puts $out [format $PDB_MASK "ATOM" $ifinta $myName $myResName $myResid $myX $myY $myZ $myOccup $myBeta $myName]
	    }
	}

	# get and write connectivity
	#BBBBB
	for {set i 1} {$i <= $NCL} {incr i} { 
	    set numbonds($i) "0"
	}
	set NWM [expr "5*$NCL"]
	for {set i 1} {$i <= $NWM} {incr i} { 
	    set nmb($i,n) "0"
	    set nmb($i,1) "0"
	    set nmb($i,2) "0"
	    set nmb($i,3) "0"
	    set nmb($i,4) "0"
	}
	set fp [open "connectivity" r]
	set file_data [read $fp] 
	set dado [split $file_data "\n"]
	foreach line $dado {
	    set j "0"
	    foreach id $line {
		set vara($j) $id
		set j [expr "$j+1"]
	    }
	    if {$TASK!="MICROSTATES"} {
		set F $CLUSTERS($vara(0),fes)
		if {$F<1000.} {
		    set myres0 "MIC"
		} else {
		    set myres0 "BAD"
		}
		set F $CLUSTERS($vara(1),fes)
		if {$F<1000.} {
		    set myres1 "MIC"
		} else {
		    set myres1 "BAD"
		}
	    } else {
		set myres0 "MIC"
		set myres1 "MIC"
	    }
	    # Just generate the bonding if all are MIC
	    if {$myres0=="MIC"} {
		if {$myres1=="MIC"} {
		    set a1 [ expr "$numbonds($vara(0)) % 5" ]
		    set a2 [ expr "$numbonds($vara(1)) % 5" ]
		    set aa1 [ expr "$a1*$NCL+$vara(0)" ]
		    set aa2 [ expr "$a2*$NCL+$vara(1)" ]
		    #
		    set nmb($aa1,n) [expr "$nmb($aa1,n)+1"]
		    set nmb($aa2,n) [expr "$nmb($aa2,n)+1"]
		    set nmb($aa1,$nmb($aa1,n)) $aa2
		    set nmb($aa2,$nmb($aa2,n)) $aa1
		    if {$nmb($aa1,n)=="4"} {
			puts $out [format $CON_MASK "CONECT" $aa1 $nmb($aa1,1) $nmb($aa1,2) $nmb($aa1,3) $nmb($aa1,4)]
			set nmb($aa1,n) "0"
		    }
		    if {$nmb($aa2,n)=="4"} {
			puts $out [format $CON_MASK "CONECT" $aa2 $nmb($aa2,1) $nmb($aa2,2) $nmb($aa2,3) $nmb($aa2,4)]
			set nmb($aa2,n) "0"
		    }
		    set numbonds($vara(0)) [expr "$numbonds($vara(0))+1"]
		    set numbonds($vara(1)) [expr "$numbonds($vara(0))+1"]
		}
	    }
	}
	set C1_MASK "%6s%5i%5i "
	set C2_MASK "%6s%5i%5i%5i "
	set C3_MASK "%6s%5i%5i%5i%5i "

	for {set i 1} {$i <= $NWM} {incr i} { 
	    if {$nmb($i,n)=="1"} {
		puts $out [format $C1_MASK "CONECT" $i $nmb($i,1) ]
	    }
	    if {$nmb($i,n)=="2"} {
		puts $out [format $C2_MASK "CONECT" $i $nmb($i,1) $nmb($i,2) ]
	    }
	    if {$nmb($i,n)=="3"} {
		puts $out [format $C3_MASK "CONECT" $i $nmb($i,1) $nmb($i,2) $nmb($i,3) ]
	    }
	}


	close $out

	puts ""
	puts "$NCL CLUSTERS, FES and BASINS written to file $CLUSTERS_and_FES_and_BASINS_PDBfile"
	puts "  X coordinate corresponds to CV number: $myVARS($Xvar,number)"
	if {$Yvar=="F"} {
	    puts "  Y coordinate corresponds to FREE ENERGY"
	} else {
	    puts "  Y coordinate corresponds to CV number: $myVARS($Yvar,number)"
	}
	if {$Zvar=="F"} {
	    puts "  Z coordinate corresponds to FREE ENERGY"
	} elseif {$Zvar=="zero"} {
	    puts "  Z coordinate set to ZERO (2D plot)"
	} else {
	    puts "  Z coordinate corresponds to CV number: $myVARS($Zvar,number)"
	}
	puts ""
    }

    #
    # loads the PDB file containing a representation of the
    # microstates and basins subdivision of the system. Called from show_basins
    proc load_fes {{CLUSTERS_and_FES_and_BASINS_PDBfile "CLUSTERS.FES.BASINS.pdb"} {TASK "BASINS"} {FES_RANGE 100}} {
	variable FES_molid
	variable TRAJ_molid
	variable START_FOLDER
	variable NCL
	cd $START_FOLDER

	if {[info exists FES_molid]} {
	    mol delete $FES_molid
	    unset FES_molid
	}

	mol new $CLUSTERS_and_FES_and_BASINS_PDBfile type pdb first 0 last -1 step 1  
	set FES_molid [molinfo top get id]

	colour_fes $FES_molid $TASK $FES_RANGE

	mol top $TRAJ_molid 

	if {[info exists pickingON]} {
	    set hola "adeu"
	} else {
	    stop_picking_clusters
	    start_picking_clusters
	}

	vmd_draw_pbcbox $FES_molid -style tubes -width 1  -color gray
    }

    # represent the microstates and basins subdivision as spheres
    # colored according to the basin they belong to. called from load_fes
    proc colour_fes {{mol_id -1} {TASK "BASINS"} {FES_RANGE 100}} {
	variable NB
	variable FES_MIN
	variable FES_MAX
	variable reprFES_MAX
	variable FES_molid
	variable nREPSinFES
	variable NCL

	if {$mol_id==-1} {
	    set mol_id $FES_molid
	}
	
	set myMIN 0
	if {[info exists reprFES_MAX]} {
	    set myMAX [expr "($reprFES_MAX-$FES_MIN)/($FES_MAX-$FES_MIN)"]
	} else {
	    set myMAX 1
	}
	clean_reps $mol_id
	# ???
	set myMIN 0
	if {[info exists reprFES_MAX]} {
	    set myMAX [expr "($reprFES_MAX-$FES_MIN)/($FES_MAX-$FES_MIN)"]
	} else {
	    set myMAX 1
	}    
	# ???
	if {$TASK=="MICROSTATES"} {

	    mol representation "CPK 3.0 1.0"
	    mol selection "all"
	    mol color "ColorID 6"
	    mol material "Opaque"
	    set rep_id [molinfo $mol_id get numreps]
	    mol addrep $mol_id
	    mol scaleminmax $mol_id $rep_id $myMIN $myMAX
	}

	if {$TASK=="FES"} {
	    set mF4 [expr "($FES_RANGE-$FES_MIN)"]
	    set mF3 [expr "($FES_RANGE-$FES_MIN)*0.75"]
	    set mF2 [expr "($FES_RANGE-$FES_MIN)*0.5"]
	    set mF1 [expr "($FES_RANGE-$FES_MIN)*0.25"]

	    mol representation "CPK 11.0 0.55"
	    mol selection "occupancy == 0"
	    mol color "Beta"
	    mol color "ColorID 1"
	    mol material "Opaque"
	    set rep_id [molinfo $mol_id get numreps]
	    mol addrep $mol_id
	    mol scaleminmax $mol_id $rep_id $myMIN $myMAX

	    mol representation "CPK 8.0 2."
	    mol selection "occupancy < $mF1"
	    mol color "Beta"
	    mol color "ColorID 1"
	    mol material "Opaque"
	    set rep_id [molinfo $mol_id get numreps]
	    mol addrep $mol_id
	    mol scaleminmax $mol_id $rep_id $myMIN $myMAX

	    mol representation "CPK 5.5 1.5"
	    mol selection "occupancy < $mF2"
	    mol color "Beta"
	    mol color "ColorID 3"
	    mol material "Opaque"
	    set rep_id [molinfo $mol_id get numreps]
	    mol addrep $mol_id
	    mol scaleminmax $mol_id $rep_id $myMIN $myMAX

	    mol representation "CPK 3. 1.3"
	    mol selection "occupancy < $mF3"
	    mol color "Beta"
	    mol color "ColorID 6"
	    mol material "Opaque"
	    set rep_id [molinfo $mol_id get numreps]
	    mol addrep $mol_id
	    mol scaleminmax $mol_id $rep_id $myMIN $myMAX

	    mol representation "CPK 2. 1."
	    mol selection "occupancy < $mF4"
	    mol color "Beta"
	    mol color "ColorID 0"
	    mol material "Opaque"
	    set rep_id [molinfo $mol_id get numreps]
	    mol addrep $mol_id
	    mol scaleminmax $mol_id $rep_id $myMIN $myMAX
	}

	if {$TASK=="BASINS"} {
	    #BBBBB
	    #       mol representation "bonds 0.1"
	    #       mol selection "resname MIC"
	    #       mol color "ColorID 0"
	    #       mol material "Opaque"
	    #       set rep_id [molinfo $mol_id get numreps]
	    #       mol addrep $mol_id
	    #       mol scaleminmax $mol_id $rep_id $myMIN $myMAX
	    set mF4 [expr "($FES_RANGE-$FES_MIN)"]

	    for {set b 0} {$b<$NB} {incr b} {
		set r [expr "$b + 1"]
		set c [expr "$b + 2"]

		mol representation "CPK 5.0 1.5"
		mol selection "(resid $r) and (occupancy < $mF4)"
		mol color "Beta"
		mol color "ColorID $c"
		mol material "Opaque"
		set rep_id [molinfo $mol_id get numreps]
		mol addrep $mol_id
		mol scaleminmax $mol_id $rep_id $myMIN $myMAX

		# attractors are represented as big spheres
		mol representation "VDW 1.5"
		mol selection "(resid $r and name C) and (occupancy < $mF4)"
		mol color "Occupancy"
		mol color "ColorID $c"
		mol material "Opaque"
		set rep_id [molinfo $mol_id get numreps]
		mol addrep $mol_id
		mol scaleminmax $mol_id $rep_id $myMIN $myMAX

	    }
	}
	set nREPSinFES [molinfo $mol_id get numreps]
	return ""
    }

    #
    # show molecular structure on VMD display
    proc show_system {} {
	variable FES_molid
	variable TRAJ_molid

	mol top $TRAJ_molid
	mol on  $TRAJ_molid
	mol off $FES_molid
	display resetview
    }

    #
    # show microstates representation on VMD display (spheres)
    proc show_fes {} {
	variable FES_molid
	variable TRAJ_molid

	mol top $FES_molid
	mol on  $FES_molid
	mol off $TRAJ_molid
	display resetview
	mol top $TRAJ_molid
    }

    #
    # show both molecular structure and microstates representation (spheres) on VMD display
    proc show_both {} {
	variable FES_molid
	variable TRAJ_molid

	mol top $TRAJ_molid
	mol on  $TRAJ_molid
	mol on  $FES_molid
	display resetview
    }

    #
    # delete all representations of the molecular structure from VMD
    proc clean_reps {{mol_id top}} {
	# clean all representations on $mol_id
	set n_reps [molinfo $mol_id get numreps]
	for {set i 0} {$i < $n_reps} {incr i} {
	    mol delrep 0 $mol_id
	}
    }

    #
    # write structures of selected microstate into a external PDB file
    proc write_clusters_frames_pdb {{cluster_id 1}} {
	variable CLUSTERS
	variable COLVAR
	variable myVARS
	variable CLUSTER_SPACE_dim
	variable TRAJ_molid
	
	#XX
        set ctrm [file exists STRUCTURES]
        if { $ctrm == 0 } {  file mkdir STRUCTURES }
	file delete -force STRUCTURES/tmp.pdb

	set sel [atomselect $TRAJ_molid "all"]
	foreach c $cluster_id {
	    file delete -force STRUCTURES/CLUSTER_${c}.pdb
	    puts ""
	    puts "writing [llength $CLUSTERS($c,frames)] structures of microstate $c into STRUCTURES/CLUSTER_${c}.pdb"
	    #    puts -nonewline "  dumping frame"
	    foreach f $CLUSTERS($c,frames) {
		#      puts -nonewline " $f"
		$sel frame $f
		$sel writepdb "STRUCTURES/tmp.pdb"
		appendfile STRUCTURES/tmp.pdb  STRUCTURES/CLUSTER_${c}.pdb
		file delete -force STRUCTURES/tmp.pdb
	    }
	    puts "done."
	    puts ""
	}
	$sel delete
	file delete -force STRUCTURES/tmp.pdb
	puts ""
    }

    # read first file, append to second, created if does not
    # exist. System-independent replacement of "exec cat fin >> fout"
    proc appendfile {fnin fnout} {
	set fin [open $fnin r]
	set fout [open $fnout a]
	puts -nonewline $fout [read $fin]
	close $fout
	close $fin
    }

    # align structures inside a cluster to frame 0
    proc align_cluster {{cluster_id 1} {salign "noh"}} {
	variable CLUSTERS
	variable COLVAR
	variable myVARS
	variable CLUSTER_SPACE_dim
	variable FRAMES
	variable TRAJ_molid

	puts " "
	puts "aligning structures inside cluster $cluster_id on selection $salign"

	set sel0 [atomselect $TRAJ_molid "$salign"]
	set sel1 [atomselect $TRAJ_molid "$salign"]
	set sel2 [atomselect $TRAJ_molid "$salign"]
	set selmv [atomselect $TRAJ_molid "all"]
	set selmv2 [atomselect $TRAJ_molid "all"]
	$sel0 frame 0
	foreach c $cluster_id {
	    set if1 [lindex $CLUSTERS($c,frames) 0]
	    $sel1 frame $if1 
	    $selmv2 frame $if1 
	    # fit first frame of the cluster to first frame of the trajectory
	    $selmv2 move [measure fit $sel1 $sel0]
	    # fit all frames of the cluster to first frame of the cluster
	    foreach f $CLUSTERS($c,frames) {
		$sel2 frame $f
		$selmv frame $f
		$selmv move [measure fit $sel2 $sel1]
	    }
	}
	puts ""
	$sel0 delete
	$sel1 delete
	$sel2 delete
	$selmv delete
    }

    #
    # writes the value of the CV for the selected microstate.
    # if structure==1 plots the  representation of its structures on the graphic window
    proc select_clusters {{cluster_id 1} {structure 0} {details 0}} {
	variable CLUSTERS
	variable COLVAR
	variable myVARS
	variable CLUSTER_SPACE_dim
	variable FRAMES
	variable TRAJ_molid
	variable MAX_FRAMES_PER_CLUSTER
	variable salign
	
	empty_repr $TRAJ_molid
	align_cluster $cluster_id $salign

	foreach c $cluster_id {
	    puts "CLUSTER $c now selected."
	    puts "  center CV coordinates: $CLUSTERS($c,center)"  
	    puts "  # elements:    $CLUSTERS($c,size)"
	    puts "  # structures:  [llength $CLUSTERS($c,frames)]"
	    puts "  Free Energy:   $CLUSTERS($c,fes)"
	    if {[info exists CLUSTERS($c,basin)]} {
		puts "  BASIN:         $CLUSTERS($c,basin)"
	    }
	    set count 0
	    foreach f $CLUSTERS($c,frames) {
		incr count
		set show 0
		if {$structure==0} {
		    if {$count<=$MAX_FRAMES_PER_CLUSTER} {
			set show 1  
		    }  
		} else {
		    foreach str $structure {
			if {$count==$str} {set show 1}
		    }
		} 
		if {$show==1} {
		    set frame_id $f
		    plot_repr $TRAJ_molid $frame_id
		    if {$details>0} {
			set ci $FRAMES($f,colvar_index)
			set cf $FRAMES($f,colvar_frame)
			puts -nonewline "    debug: frame $f coords / walker / colvar frame / bias : \t"
			for {set d 0} {$d < $CLUSTER_SPACE_dim} {incr d} {
			    puts -nonewline "$COLVAR($ci,$cf,$myVARS($d,number)) "
			}
			puts -nonewline "$ci "
			puts -nonewline "$cf "
			puts -nonewline "$COLVAR($ci,$cf,bias) "
			puts ""
		    }
		}
	    }
	}
	if {[info exists frame_id]} {
	    animate goto $frame_id
	}
    }

    #
    # writes the value of the CV for the selected basin and represents its structures on the graphic window
    proc select_basin {{basin_id 1}} {
	variable BASINS
	variable CLUSTERS

	select_clusters $BASINS($basin_id,clusters_list)

	puts ""
	puts "BASIN $basin_id now selected."
	puts "  cluster center: $BASINS($basin_id,cluster_center)"  
	puts "  center CV coordinates: $CLUSTERS($BASINS($basin_id,cluster_center),center)"  
	puts "  # clusters:       $BASINS($basin_id,size)"
	set strs 0
	set els 0
	foreach cc $BASINS($basin_id,clusters_list) {
	    set str [llength $CLUSTERS($cc,frames)]
	    set el $CLUSTERS($cc,size)
	    set strs [expr "$strs + $str"]
	    set els [expr "$els + $el"]
	}
	puts "  # elements:       $els"
	puts "  # structures:     $strs"
	puts "  Av. Free Energy:  $BASINS($basin_id,avrg_fes)"
    }

    #
    # plots the  representation of the structures in  frame_id on the graphic window
    proc plot_repr {mol_id frame_id} {
	variable FES_MIN
	variable FES_MAX
	variable reprFES_MAX

	set myMIN $FES_MIN
	if {[info exists reprFES_MAX]} {
	    set myMAX $reprFES_MAX
	} else {
	    set myMAX $FES_MAX
	}

	set n_reps [molinfo $mol_id get numreps]
	for {set rep_id 0} {$rep_id < $n_reps} {incr rep_id} {
	    set frames_list [mol drawframes $mol_id $rep_id]
	    lappend frames_list $frame_id
	    mol drawframes $mol_id $rep_id "$frames_list"
	    mol scaleminmax $mol_id $rep_id $myMIN $myMAX
	}
    }

    proc empty_repr {mol_id} {
	set n_reps [molinfo $mol_id get numreps]
	for {set rep_id 0} {$rep_id < $n_reps} {incr rep_id} {
	    mol drawframes $mol_id $rep_id ""
	    mol scaleminmax $mol_id $rep_id 0.00 0.00
	}
    }

    # 
    # reads the FES file, containing the free energy; sets the variable CLUSTERS($cluster_id,fes)
    proc read_fes {{FES_file "WHAM_RUN/FES"}} {
	variable CLUSTERS
	variable clusters_LIST
	variable ncLIST
	variable FES_MIN
	variable FES_MAX
	variable NCL
	variable START_FOLDER
	cd $START_FOLDER

	puts "resetting all clusters fes to 1000"
	for {set i 1} {$i<=$NCL} {incr i} {
	    set CLUSTERS($i,fes) 1000
	}

	set FES_MIN 999999999
	set FES_MAX -999999999

	puts "reading fes from $FES_file"
	set i 0
	set in [open $FES_file r]
	while {[gets $in line] != -1} {
	    set cluster_id [lindex $line 0]
	    set fes [lindex $line 1]
	    set CLUSTERS($cluster_id,fes) $fes
	    incr i
	    set clusters_LIST($i) $cluster_id
	    if {$fes < $FES_MIN} {set FES_MIN $fes} 
	    if {$fes < 998. && $fes > $FES_MAX} {set FES_MAX $fes} 
	}
	set ncLIST $i

	puts "$i cluster fes read"
    }

    #
    # assigns microstate free energy to each structure in the varaible "user" (for coloring the structures according to the free energy)
    proc put_fes_on_user {} {
	variable CLUSTERS
	variable FRAMES
	variable TRAJ_molid
	variable TOTAL_TRAJ_FRAMES

	puts "assigning FES to each frame of the TRAJECTORY"
	set sel [atomselect $TRAJ_molid "all"]
	for {set f 0} {$f < $TOTAL_TRAJ_FRAMES} {incr f} {
	    set fes $CLUSTERS($FRAMES($f,cluster),fes)
	    $sel frame $f
	    $sel set user $fes
	}

	$sel delete
	puts "done"
	puts ""
    }

    #
    # reads a file containing a list of microstates. Sets the variable clusters_LIST($i) 
    proc read_clusters_list {clusters_LIST_file} {
	variable clusters_LIST
	variable ncLIST
	variable START_FOLDER
	cd $START_FOLDER

	set i 0
	set in [open $clusters_LIST_file r]
	while {[gets $in line] != -1} {
	    foreach cluster_id $line {
		incr i
		set clusters_LIST($i) $cluster_id
	    }
	}
	set ncLIST $i
	puts "$i cluster ids read"
    }

    proc clean_fes {} {
	variable nREPSinFES
	variable FES_molid

	set mol_id $FES_molid
	set n_reps [molinfo $mol_id get numreps]
	for {set i $n_reps} {$i >= $nREPSinFES} {set i [expr "$i - 1"]} {
	    mol delrep $i $mol_id
	}
	label delete Atoms all
    }

    #
    # creates a large sphere on the microstate that has been selected
    proc highlight_fes {cluster_id} {
	variable FES_molid
	variable FES_MIN
	variable FES_MAX
	variable reprFES_MAX

	set mol_id $FES_molid

	set myMIN 0
	if {[info exists reprFES_MAX]} {
	    set myMAX [expr "($reprFES_MAX-$FES_MIN)/($FES_MAX-$FES_MIN)"]
	} else {
	    set myMAX 1
	}

	mol representation "VDW 2.0"
	mol color "Beta"
	mol material "Opaque"
	mol selection "serial $cluster_id"
	set rep_id [molinfo $mol_id get numreps]
	mol addrep $mol_id
	mol selupdate $rep_id $mol_id "on"
	mol scaleminmax $mol_id $rep_id $myMIN $myMAX
    }

    #
    # enables picking microstates with the mouse. called from load_fes
    proc start_picking_clusters {} {
	variable pickingON
	trace add variable ::vmd_pick_event write ::metagui3::show_cluster_from_pick
	set pickingON 1
	puts "picking clusters from free energy is now ON"
    }

    #
    # disables picking microstates with the mouse. called from load_fes
    proc stop_picking_clusters {} {
	variable pickingON
	trace vdelete ::vmd_pick_event w ::metagui3::show_cluster_from_pick
	if {[info exists pickingON]} {
	    unset pickingON
	    puts "picking clusters from free energy is now OFF"
	}
    }

    #
    # driver function for picking microstates with the mouse. Calls
    # select_clusters for the picked microstate and highlights the
    # selected microstate
    proc show_cluster_from_pick {args} {
	global vmd_pick_atom vmd_pick_mol
	variable TRAJ_molid
	variable BASINS
	variable FES_molid

	# use the picked atom's index and molecule id
	puts "In show_cluster_from_pick, $args"

	if {$vmd_pick_mol==$FES_molid} {
	    set sel [atomselect $vmd_pick_mol "index $vmd_pick_atom"]
	    set basin_id [$sel get resid]
	    set cluster_id [$sel get serial]
	    #set cluster_id $BASINS($basin_id,cluster_center)
	    
	    select_clusters $cluster_id
	    
	    clean_fes
	    highlight_fes $cluster_id 

	    show_system
	}
    }

    #
    # produces multiplot graphs with the history-dependent potential used in the WHAM analysis (one-dimensional bias)
    proc show_wham_graphs {{hill_id 1} {WHAM_FOLDER "WHAM_RUN"}} {
	variable CV_INFO_GRID
	variable START_FOLDER
	variable TRAJ_INFO
	variable plothandle
	variable ploth
	variable winf
	cd $START_FOLDER
	set i [expr "($hill_id + 1)"]
	set bn [file tail $TRAJ_INFO([hills2traj $hill_id],hills_file)]
	set INPUT_file [format "VG_%s" $bn]
	set CVNAME $CV_INFO_GRID($i,type)
	cd $WHAM_FOLDER
	if { [winfo exists $winf] == 1 } {
	    $ploth quit
	}
	set in [open $INPUT_file r]
	set x ""
	set xav ""
	set xa ""
	set xc ""
	set firstp ""
	set secondp ""
	set avg ""
	set avgwr ""
	set conn ""
	while {[gets $in line] != -1} {
	    if {[llength $line] > 0} {
		set c1 [lindex $line 0]
		set c2 [lindex $line 1]
		set c3 [lindex $line 2]
		set c4 [lindex $line 3]
		set c5 [lindex $line 4]
		set c6 [lindex $line 5]
		set c2 [expr "($c2 * -1)"]
		set c3 [expr "($c3 * -1)"]
		set c4 [expr "($c4 * -1)"]
		set c5 [expr "($c5 * -1)"]
		set c6 [expr "($c6 * -1)"]
		lappend x $c1
		lappend firstp $c2
		lappend secondp $c3
		if {$c4 < 998} {
		    lappend xav $c1
		    lappend avg $c4
		}
		if {$c5 < 998} {
		    lappend xa $c1
		    lappend avgwr $c5
		}
		if {$c6 < 998} {
		    lappend xc $c1
		    lappend conn $c6
		}
	    }
	}
	set title $INPUT_file 
	set ploth [multiplot -x $x -y $firstp -xsize 500 -ysize 375 -title $title -lines -linewidth 1 -linecolor red -legend "first part"  -ylabel "F" -xlabel $CVNAME -plot]
	$ploth add $x $secondp -lines -linewidth 1  -linecolor blue -legend "second part" -plot
	$ploth add $xav $avg -lines -linewidth 2  -linecolor black -legend "average" -plot
	$ploth add $xa $avgwr -nolines -marker circle -radius 6 -fillcolor black -legend "within range" -plot
	$ploth add $xc $conn -nolines -marker circle -radius 4 -fillcolor green -legend "connected" -plot
	$ploth replot;
	set winf [$ploth getpath]
	cd ..
    }

    # produces multiplot graph with time scales for the kinetic basins analysis
    proc show_TAU {{BASINS_FOLDER "BASINS_RUN"}} {
	variable START_FOLDER
	variable plothandle
	variable ploth
	variable winf
	cd $START_FOLDER
	set INPUT_file "TAU"
	cd $BASINS_FOLDER
	if { [winfo exists $winf] == 1 } {
	    $ploth quit
	}
	set in [open $INPUT_file r]
	set x ""
	set firstp ""
	set secondp ""
	set avg ""
	set avgwr ""
	set conn ""
	while {[gets $in line] != -1} {
	    if {[llength $line] > 0} {
		set c1 [lindex $line 0]
		set c2 [lindex $line 1]
		lappend x $c1
		lappend firstp $c2
	    }
	}
	set title $INPUT_file 
	set ploth [multiplot -x $x -y $firstp -xsize 500 -ysize 375 -title $title -nolines -marker circle -radius 6 -fillcolor red -legend ""  -ylabel "TAU" -xlabel "N" -plot]
	$ploth replot;
	set winf [$ploth getpath]
	cd ..
    }

    # produces multiplot graph with decision graph for smart clustering
    proc show_DEC {{SAF_FOLDER "SAF_RUN"}} {
	variable START_FOLDER
	variable plothandle
	variable ploth
	variable winf
	cd $START_FOLDER
	set INPUT_file "dec.dat"
	cd $SAF_FOLDER
	if { [winfo exists $winf] == 1 } {
	    $ploth quit
	}
	set in [open $INPUT_file r]
	set x "";		# Rho vector
	set firstp "";		# Delta vector
	set secondp ""
	set avg ""
	set avgwr ""
	set conn ""
	while {[gets $in line] != -1} {
	    if {[llength $line] > 0} {
		lappend x [lindex $line 1]
		lappend firstp [lindex $line 2]
	    }
	}
	set xsort [lsort -real $x]
	set xmin [lindex $xsort 0]
	set xmax [lindex $xsort end]
	set dis [expr {($xmax - $xmin) /10.} ]
	set dis [format "%.2f" $dis]
	puts "INFO: rho limits are xmin=$xmin xmax=$xmax dis=$dis"
	set title $INPUT_file 
	set ploth [multiplot -x $x -y $firstp -title $title -xsize 500 -ysize 375 \
		       -nolines  -marker circle -radius 6 \
		       -fillcolor red -legend ""  -ylabel $::metagui3::u8_delta \
		       -xlabel  $::metagui3::u8_rho ]
	# TONI for some reason -xmajortics $dis causes a "too many
	# nested evaluations" error
	$ploth replot;
	set winf [$ploth getpath]
	# TONI the following assumes that "find withtag point" returns
	# the points in the same order as they were inserted, and adds
	# bindings for left and right mouse clicks.
	foreach id [$winf.f.cf find withtag point] rho $x delta $firstp {
	    $winf.f.cf bind $id <Button-1> "set ::metagui3::RHO_CUT $rho"
	    $winf.f.cf bind $id <Button-3> "set ::metagui3::DELTA_CUT $delta"
	    $winf.f.cf bind $id <Control-Button-1> "set ::metagui3::DELTA_CUT $delta"	
#	    $winf.f.cf bind $id <Button-1> "::metagui3::show_DEC_cut rho $rho $ploth"
#	    $winf.f.cf bind $id <Button-3> "::metagui3::show_DEC_cut delta $delta $ploth" 
#	    $winf.f.cf bind $id <Control-Button-1>  "::metagui3::show_DEC_cut delta $delta $ploth" 
	}
	cd ..
    }

    # Procedures to set and show rho/delta cuts. For some reason the
    # bindings are lost after the first click!?!
    proc show_DEC_cut {axis value ploth} {
	variable RHO_CUT
	variable DELTA_CUT

	switch $axis {
	    rho {
		set RHO_CUT $value
		$ploth configure -vline [list $value]
	    }
	    delta {
		set DELTA_CUT $value
		$ploth configure -hline [list $value]
	    }
	}
	$ploth replot
    }


    # Return the index in TRAJ_INFO (formerly COLVAR_INFO)
    # corresponding t the N-th biased trajectory (both 0-based). This
    # is done to avoid having separate and duplicate HILLS and
    # HILLS_INFO arrays.
    proc hills2traj {N} {
	variable TRAJ_INFO
	variable nTRAJECTORIES
	set nbias 0
	for {set i 0} {$i < $nTRAJECTORIES} {incr i} {
	    if {$TRAJ_INFO($i,bias)==1} {
		if {$N==$nbias} {
		    return $i
		}
		incr nbias
	    }
	}
	error "Internal error at hills2traj: hills $N not found"
    }

    # Returns the number of active hills
    proc nHILLS {} {
	variable TRAJ_INFO
	variable nTRAJECTORIES
	set nbias 0
	for {set i 0} {$i < $nTRAJECTORIES} {incr i} {
	    if {$TRAJ_INFO($i,bias)==1} {
		incr nbias
	    }
	}
	return $nbias
    }

    # changes the list of active CVs (CV_ACTIVE) when one clicks on a
    # "use" button in the CV section. TONI making this explicit is
    # redundant and error prone.  The information is already in CV
    # GRID, duplicated in nCV_ACTIVE, and potentially outdated.
    proc compute_active {} {
	variable CV_INFO_GRID
	variable nVARIABLES
	variable nCV_ACTIVE
	variable CV_ACTIVE
	variable out2

	set nCV_ACTIVE 0
	array unset CV_ACTIVE

	for { set i  1} {$i <= $nVARIABLES} { incr i } {
	    set CV_INFO_GRID($i,useplot) "No"
	    if {$CV_INFO_GRID($i,use)==1} {
		set CV_ACTIVE($nCV_ACTIVE) $i
		incr nCV_ACTIVE 
	    }
	}
    }


    
}



