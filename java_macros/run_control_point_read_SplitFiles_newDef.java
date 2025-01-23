// Simcenter STAR-CCM+ macro: run_control_point_read.java
// Written by Michael Barbour
// Macro manages a simulation with control point defined motion. At a specific time interval, a new contorl points file is read in and the displacement field is updated

package macro;

import java.util.*;
import java.io.*;
import java.nio.file.*;

import star.common.*;
import star.morpher.*;

public class run_control_point_read_SplitFiles_newDef extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    int n_periods = 150;
    int period_count = 120;  //start
    Double period_length = 0.01; //spacing in each cp file

    File f = null;
    String[] paths;

    f = new File("D:\\Barbour\\OneDrive - UW\\TracheaMalacia\\CFD\\Motion_test\\Small_start\\new_cp_test_0001\\tmp\\");
    paths = f.list();


    Simulation simulation_0 = 
      getActiveSimulation();

    PhysicalTimeStoppingCriterion stop_time = 
      ((PhysicalTimeStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Physical Time"));


    while(period_count <= (n_periods-1)){


      // Get the correct file table 
      String table_path = paths[period_count];
      String[] parts = table_path.split("\\.(?=[^\\.]+$)");
      String table_name = parts[0];

      simulation_0.println("Loading CP Table: " + table_name);

      // Update PointSets displacement field with correct table
      PointSet pointSet_0 = 
        ((PointSet) simulation_0.get(PointSetManager.class).getObject("Table Point Set"));

      // TablePointGenerator tablePointGenerator_0 = 
      //   ((TablePointGenerator) pointSet_0.getPointGenerator());

      FileTable fileTable_100 = 
         ((FileTable) simulation_0.getTableManager().getTable(table_name));

      // tablePointGenerator_0.setTable(fileTable_100);

      // tablePointGenerator_0.setX0Data("X");

      // tablePointGenerator_0.setY0Data("Y");

      // tablePointGenerator_0.setZ0Data("Z");

      // tablePointGenerator_0.regeneratePointSet();

      IncrementalDisplacementProfile incrementalDisplacementProfile_0 = 
        pointSet_0.getValues().get(IncrementalDisplacementProfile.class);

      incrementalDisplacementProfile_0.getMethod(TimeXyzTabularVectorProfileMethod.class).setTable(fileTable_100);


      //Set max physical time equal to period length and run until then
      period_count = period_count + 1;
      stop_time.setMaximumTime(period_length * period_count);
      simulation_0.getSimulationIterator().run();

      Double time = getActiveSimulation().getSolution().getPhysicalTime();
      simulation_0.println("Simulation Time: " + time);



    }
  }
}
