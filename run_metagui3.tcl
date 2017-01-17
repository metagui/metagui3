# Convenience script to start up metagui

# Load TK shell on sourcing
menu tkcon on

# Should be done by the user, but let's do it again in case they
# didn't, to be sure. I know no reliable way to say "the dir this
# script is in". This won't work with e.g. "vmd -e
# ../run_metagui3.tcl". Workaround: TCLLIBPATH=... vmd -e ...
lappend auto_path [file dirname [info script]]

# Loads metagui functions
package require metagui3 

# Launch Graphical User Interface immediately after sourcing
metagui3::metagui

# Load metagui.in if exists
if {[file exists "metagui.in"]} {
    metagui3::read_input "metagui.in"
}

# If not done already - perhaps this should go in pkgIndex.tcl
vmd_install_extension metagui3 metagui3::metagui "Analysis/Metagui 2.0"

