// Simcenter STAR-CCM+ macro: read_control_points.java
// Written by Michael Barbour, Simcenter STAR-CCM+ 18.06.006
// Reads in all control point files located in a specified directory
//

package macro;

import java.util.*;
import java.io.*;
import java.nio.file.*;

import star.common.*;
import star.base.neo.*;


public class read_control_points extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    Simulation simulation_0 = 
      getActiveSimulation();

    File f = null;
    String[] paths;
    int numFiles = 0;


    f = new File("D:\\Barbour\\OneDrive - UW\\TracheaMalacia\\CFD\\Motion_test\\Small_start\\Split_CPs_fullTime\\");
    paths = f.list();


    for(String path:paths){

      simulation_0.println(path);

      FileTable fileTable_100 = (FileTable) simulation_0.getTableManager().createFromFile(resolvePath(f + "/" + path), null);
      numFiles++;
    }

    simulation_0.println("Loaded " + numFiles + " Control Point Files");
    
  }
}
