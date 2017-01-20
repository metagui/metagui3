(Alanine)3 Demo for Metagui
===========================

Usage - assuming that METAGUI3 is installed properly:

1. Start the plugin.

        cd DEMO_ALA3
        TCLLIBPATH=.. vmd -e ../run_metagui3.tcl


2. **Load the trajectories**. click "Load all"

3. **Cluster into microstates**. Go in the *Analysis/Microstates* tab
   and click on "Find microstates".

4. **Compute the free energies of each cluster.** Go in the
   *Analysis/Thermodynamics* tab and click on "Compute free energies".

5. **Inspect microstates.** The *Visualization/Microstates* table
   should be now filled and can be visualized.

6. **Group clusters into kinetic basins.** Go in the
   *Analysis/Clustering* tab and:
     * Select "Show decision plane"
     * Identify a cluster of separated points on the top-right corner.
	 * Left-click and right-click on the isolated point at the bottom-left.
	 * Select "Find kinetic basins"

7. **Inspect kinetic basins.** The *Visualization/Basins* list table
   should be now filled and can be inspected.
