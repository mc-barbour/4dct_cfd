// Simcenter STAR-CCM+ macro: parabolic_flow_profile_contd.java
// Written by Simcenter STAR-CCM+ 16.04.012
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;

public class parabolic_flow_profile_contd extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    Simulation simulation_0 = 
      getActiveSimulation();

    UserFieldFunction userFieldFunction_1 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_1.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SCALAR);

    userFieldFunction_1.setPresentationName("trachea_radius");

    userFieldFunction_1.setFunctionName("trachea_radius");

    userFieldFunction_1.setDimensions(Dimensions.Builder().length(1).build());

    userFieldFunction_1.setDefinition("0.0022977");


    UserFieldFunction userFieldFunction_2 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_2.getTypeOption().setSelected(FieldFunctionTypeOption.Type.POSITION);

    userFieldFunction_2.setFunctionName("trachea_center");

    userFieldFunction_2.setPresentationName("trachea_center");

    userFieldFunction_2.setDimensions(Dimensions.Builder().length(1).build());

    userFieldFunction_2.setDefinition("[-0.006836975226178765, -0.11067786812782288, 0.0031251306645572186]");

    UserFieldFunction userFieldFunction_3 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_3.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SCALAR);

    userFieldFunction_3.setPresentationName("radial_position");

    userFieldFunction_3.setFunctionName("radial_position");

    userFieldFunction_3.setDimensions(Dimensions.Builder().length(1).build());

    userFieldFunction_3.setDefinition("sqrt(pow($${Position}[0] - $${trachea_center}[0], 2) + pow($${Position}[1] - $${trachea_center}[1], 2) + pow($${Position}[2] - $${trachea_center}[2],2))  ");

    UserFieldFunction userFieldFunction_4 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_4.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SCALAR);

    userFieldFunction_4.setPresentationName("velocity_max");

    userFieldFunction_4.setFunctionName("velocity_max");

    userFieldFunction_4.setDimensions(Dimensions.Builder().length(1).time(-1).build());

    userFieldFunction_4.setDefinition("2. * 8.47E-05 / ${Density} / (3.14 * pow(${trachea_radius},2)) ");

    UserFieldFunction userFieldFunction_5 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_5.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SCALAR);

    userFieldFunction_5.setPresentationName("velocity_profile");

    userFieldFunction_5.setFunctionName("velocity_profile");

    userFieldFunction_5.setDimensions(Dimensions.Builder().length(1).time(-1).build());

    userFieldFunction_5.setDefinition("${velocity_max}*(1 - pow(${radial_position} / ${trachea_radius}, 2))");
  }
}
