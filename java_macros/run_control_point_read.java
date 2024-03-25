// Simcenter STAR-CCM+ macro: run_control_point_read.java
// Written by Michael Barbour
// Macro manages a simulation with control point defined motion. At a specific time interval, a new contorl points file is read in and the displacement field is updated

package macro;

import java.util.*;

import star.common.*;
import star.morpher.*;

public class run_control_point_read extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    int n_periods = 3;
    int period_count = 1;
    Double period_length = 0.4;
    String control_point_file_prefix = "cylinder_control_points_0.01_incDisp_";

    Simulation simulation_0 = 
      getActiveSimulation();

    PhysicalTimeStoppingCriterion stop_time = 
      ((PhysicalTimeStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Physical Time"));


    while(period_count <= n_periods){

      //Set max physical time equal to period length and run until then
      stop_time.setMaximumTime(period_length * period_count);
      simulation_0.getSimulationIterator().run();

      Double time = getActiveSimulation().getSolution().getPhysicalTime();
      simulation_0.println("Simulation Time: " + time);

      period_count = period_count + 1;

      if (period_count == n_periods){
        break; // Don't try and read control points file at last period
      }

      // Define and read control points file
      String control_point_filename = control_point_file_prefix + Integer.toString(period_count) + ".csv";
      simulation_0.println(control_point_filename);

      FileTable fileTable_0 = 
        ((FileTable) simulation_0.getTableManager().createFromFile(resolvePath(control_point_filename),null));

      // Update PointSets displacement field with new table
      PointSet pointSet_0 = 
      ((PointSet) simulation_0.get(PointSetManager.class).getObject("Table Point Set Run"));

      IncrementalDisplacementProfile incrementalDisplacementProfile_0 = 
      pointSet_0.getValues().get(IncrementalDisplacementProfile.class);

      incrementalDisplacementProfile_0.getMethod(TimeXyzTabularVectorProfileMethod.class).setTable(fileTable_0);



    }
  }
}
