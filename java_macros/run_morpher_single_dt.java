// Simcenter STAR-CCM+ macro: run_morpher_dingle_dt.java
// Written by Michael Barbour
// Macro manages a simulation with control point defined motion. At a specific time interval, a new contorl points file is read in and the displacement field is updated

package macro;

import java.util.*;
import java.io.*;
import java.nio.file.*;

import star.common.*;
import star.base.neo.*;
import star.morpher.*;
import star.motion.*;

public class run_morpher_single_dt extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {
    
    Double end_time = 3.0;
    Double dt = 0.01;
    int n_time_steps = (int) Math.round(end_time / dt);
    int time_step_count = 0;


    // String control_point_file_prefix = "StarControlPoinstFull_2mmSpacing_Periodic_5Cycles_Inc_interpolated_ptj05_01s_cycle";

    File f = null;
    String[] paths;

    f = new File("D:\\Barbour\\OneDrive - UW\\TracheaMalacia\\CFD\\Ama_demo_files\\cp_files_everyDT\\");
    paths = f.list();


    Simulation simulation_0 = 
      getActiveSimulation();

    simulation_0.println(n_time_steps);

    PhysicalTimeStoppingCriterion stop_time = 
      ((PhysicalTimeStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Physical Time"));
    
  
      
    // Custom comparator to compare the numeric part of the strings
    Comparator<String> numericComparator = new Comparator<String>() {
        @Override
        public int compare(String s1, String s2) {
            int num1 = extractNumber(s1);
            int num2 = extractNumber(s2);
            return Integer.compare(num1, num2);
        }

        // Method to extract the last numeric part from the string
        private int extractNumber(String s) {
            String num = s.substring(s.lastIndexOf('e') + 1, s.length() - 4); // Remove ".csv"
            return Integer.parseInt(num);
        }
    };
    
    simulation_0.println(paths[0]);
      
    Double time = getActiveSimulation().getSolution().getPhysicalTime();
    long round_time = Math.round(time/dt);
    time_step_count = (int) round_time; 
    simulation_0.println("File Index: " + time_step_count);
    
    
    
    

    while(time_step_count <= n_time_steps-1){
    
    

      String path = paths[time_step_count];
      simulation_0.println(path);
      
      //load cp excel file
      FileTable fileTable_100 = (FileTable) simulation_0.getTableManager().createFromFile(resolvePath(f + "/" + path), null);
        
      // Get the correct file table 
      String table_path = paths[time_step_count];
      String[] parts = table_path.split("\\.(?=[^\\.]+$)");
      String table_name = parts[0];

      if (time_step_count == 0){
          PointSet pointSet_0 = 
            simulation_0.get(PointSetManager.class).createTablePointSet("Table Point Set", fileTable_100, "X", "Y", "Z");
      }
      
      // Update PointSets displacement field with new table
      PointSet pointSet_1 = 
      ((PointSet) simulation_0.get(PointSetManager.class).getObject("Table Point Set"));
      
      
      TablePointGenerator tablePointGenerator_0 = 
        ((TablePointGenerator) pointSet_1.getPointGenerator());

      tablePointGenerator_0.setTable(fileTable_100);

      tablePointGenerator_0.setX0Data("X");

      tablePointGenerator_0.setY0Data("Y");

      tablePointGenerator_0.setZ0Data("Z");

      tablePointGenerator_0.regeneratePointSet();
     


      PointSetMotionSpecification pointSetMotionSpecification_0 = 
        pointSet_1.getValues().get(PointSetMotionSpecification.class);

      MorphingMotion morphingMotion_0 = 
        ((MorphingMotion) simulation_0.get(MotionManager.class).getObject("Morphing"));

      pointSetMotionSpecification_0.setMotion(morphingMotion_0);

      IncrementalDisplacementProfile incrementalDisplacementProfile_0 = 
        pointSet_1.getValues().get(IncrementalDisplacementProfile.class);

      incrementalDisplacementProfile_0.setMethod(TimeXyzTabularVectorProfileMethod.class);

      incrementalDisplacementProfile_0.getMethod(TimeXyzTabularVectorProfileMethod.class).setTable(fileTable_100);

      incrementalDisplacementProfile_0.getMethod(TimeXyzTabularVectorProfileMethod.class).setXData("X");

      incrementalDisplacementProfile_0.getMethod(TimeXyzTabularVectorProfileMethod.class).setYData("Y");

      incrementalDisplacementProfile_0.getMethod(TimeXyzTabularVectorProfileMethod.class).setZData("Z");
      /*
      FileTable to_delete = fileTable_100;
      */
      
      if (time_step_count != 0){
      
          // Get the correct file table 
          String table_path_to_delete = paths[time_step_count-1];
          String[] parts_to_delete = table_path_to_delete.split("\\.(?=[^\\.]+$)");
          String table_name_to_delete = parts_to_delete[0];
          
          simulation_0.println("Deleting Table: " + table_name_to_delete);
          
          FileTable to_delete = 
            ((FileTable) simulation_0.getTableManager().getTable(table_name_to_delete));
            
      
          simulation_0.deleteObjects(new ArrayList<>(Arrays.<ClientServerObject>asList(to_delete)));
      }
      
      //Set max physical time equal to period length and run until then
      simulation_0.println("Loading CP Table: " + table_name);
      time_step_count = time_step_count + 1;
      stop_time.setMaximumTime(dt * time_step_count);
      simulation_0.getSimulationIterator().run();
      
      time = getActiveSimulation().getSolution().getPhysicalTime();
      
      simulation_0.println("Simulation Time: " + time);
      



    }
  }
}