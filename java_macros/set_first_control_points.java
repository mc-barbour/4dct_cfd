// Simcenter STAR-CCM+ macro: control_point_set.java
// Written by Simcenter STAR-CCM+ 18.06.006
package macro;

import java.util.*;
import java.io.*;

import star.common.*;
import star.base.neo.*;
import star.morpher.*;
import star.motion.*;

public class set_first_control_points extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    // String propFileName = "D:\\Barbour\\OneDrive - UW\\RobinSequence\\Data\\CFD\\MovingMesh\\cylinder_test\\sim.properties";
    // Properties prop = new Properties();
    
    // ObjectInputStream inputStream = new ObjectInputStream(new FileInputStream(propFileName));
    // // prop.load(inputStream);
    // if (inputStream != null) {
    //   prop.load(inputStream);
    // } else {
    //   throw new FileNotFoundException("property file '" + propFileName + "' not found in the classpath");
    // }

    String control_point_file_prefix = prop.getProperty("control_point_file_prefix");

    Simulation simulation_0 = 
      getActiveSimulation();

    long start = System.nanoTime();

    FileTable fileTable_0 = 
      (FileTable) simulation_0.getTableManager().createFromFile(resolvePath("cylinder_control_points_0.01_incDisp_1.csv"), null);

    long finish = System.nanoTime();
    long timeElapsed = finish - start;
    System.out.println(timeElapsed);

    
    PointSet pointSet_0 = 
      simulation_0.get(PointSetManager.class).createTablePointSet("Table Point Set Run", fileTable_0, "X", "Y", "Z");

    TablePointGenerator tablePointGenerator_0 = 
      ((TablePointGenerator) pointSet_0.getPointGenerator());

    tablePointGenerator_0.setTable(fileTable_0);

    tablePointGenerator_0.setX0Data("X");

    tablePointGenerator_0.setY0Data("Y");

    tablePointGenerator_0.setZ0Data("Z");

    tablePointGenerator_0.regeneratePointSet();

    PointSetMotionSpecification pointSetMotionSpecification_0 = 
      pointSet_0.getValues().get(PointSetMotionSpecification.class);

    MorphingMotion morphingMotion_0 = 
      ((MorphingMotion) simulation_0.get(MotionManager.class).getObject("Morphing"));

    pointSetMotionSpecification_0.setMotion(morphingMotion_0);

    IncrementalDisplacementProfile incrementalDisplacementProfile_0 = 
      pointSet_0.getValues().get(IncrementalDisplacementProfile.class);

    incrementalDisplacementProfile_0.setMethod(TimeXyzTabularVectorProfileMethod.class);

    incrementalDisplacementProfile_0.getMethod(TimeXyzTabularVectorProfileMethod.class).setTable(fileTable_0);

    incrementalDisplacementProfile_0.getMethod(TimeXyzTabularVectorProfileMethod.class).setXData("X");

    incrementalDisplacementProfile_0.getMethod(TimeXyzTabularVectorProfileMethod.class).setYData("Y");

    incrementalDisplacementProfile_0.getMethod(TimeXyzTabularVectorProfileMethod.class).setZData("Z");
  }
}
