// Simcenter STAR-CCM+ macro: viscous_dissipation.java
// Written by Simcenter STAR-CCM+ 16.02.008
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.base.report.*;

public class viscous_dissipation extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    Simulation simulation_0 = 
      getActiveSimulation();

    UserFieldFunction userFieldFunction_1 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_1.getTypeOption().setSelected(FieldFunctionTypeOption.Type.VECTOR);

    userFieldFunction_1.setPresentationName("grad_vel_u");

    userFieldFunction_1.setFunctionName("grad_vel_u");

    userFieldFunction_1.setDefinition("grad($$Velocity[0])");

    userFieldFunction_1.setDimensions(Dimensions.Builder().time(-1).build());

    UserFieldFunction userFieldFunction_2 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_2.getTypeOption().setSelected(FieldFunctionTypeOption.Type.VECTOR);

    userFieldFunction_2.setPresentationName("grav_vel_v");

    userFieldFunction_2.setFunctionName("grad_vel_v");

    userFieldFunction_2.setDimensions(Dimensions.Builder().time(-1).build());

    userFieldFunction_2.setDefinition("grad($$Velocity[1])");

    UserFieldFunction userFieldFunction_3 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_3.getTypeOption().setSelected(FieldFunctionTypeOption.Type.VECTOR);

    userFieldFunction_3.setPresentationName("grad_vel_w");

    userFieldFunction_3.setFunctionName("grad_vel_w");

    userFieldFunction_3.setDefinition("grad($$Velocity[2])");

    UserFieldFunction userFieldFunction_4 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_4.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SYMMETRIC_TENSOR);

    userFieldFunction_4.setPresentationName("velocity_gradient_tensor");

    userFieldFunction_4.setFunctionName("velocity_gradient_tensor");

    userFieldFunction_4.setDimensions(Dimensions.Builder().time(-1).build());

    userFieldFunction_4.setDefinition("[$${grad_vel_u}[0], $${grad_vel_u}[1], $${grad_vel_u}[2]; $${grad_vel_v}[1], $${grad_vel_v}[2]; $${grad_vel_w}[2]]");

    UserFieldFunction userFieldFunction_5 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_5.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SYMMETRIC_TENSOR);

    userFieldFunction_5.setFunctionName("velocity_gradient_transpos");

    userFieldFunction_5.setFunctionName("velocity_gradient_transpose");

    userFieldFunction_5.setPresentationName("velocity_gradient_transpose");

    userFieldFunction_5.setDimensions(Dimensions.Builder().time(-1).build());

    userFieldFunction_5.setDefinition("[$${grad_vel_u}[0], $${grad_vel_v}[0], $${grad_vel_w}[0]; $${grad_vel_v}[1], $${grad_vel_w}[1]; $${grad_vel_w}[2]]");

    userFieldFunction_3.setDimensions(Dimensions.Builder().time(-1).build());

    UserFieldFunction userFieldFunction_6 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_6.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SYMMETRIC_TENSOR);

    userFieldFunction_6.setPresentationName("stress_tensor");

    userFieldFunction_6.setFunctionName("stress_tensor");

    userFieldFunction_6.setDefinition("$$$velocity_gradient_tensor + $$$velocity_gradient_transpose");

    userFieldFunction_6.setDimensions(Dimensions.Builder().time(-1).build());

    UserFieldFunction userFieldFunction_7 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_7.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SYMMETRIC_TENSOR);

    userFieldFunction_7.setPresentationName("viscous_dissipation");

    userFieldFunction_7.setFunctionName("viscous_dissipation");

    userFieldFunction_7.setDimensions(Dimensions.Builder().mass(1).length(-1).time(3).build());

    userFieldFunction_7.setDimensions(Dimensions.Builder().mass(1).length(-1).time(-3).build());

    userFieldFunction_7.setDefinition("dot($$$stress_tensor,$$$stress_tensor)*${DynamicViscosity} * 0.5");

    VolumeIntegralReport volumeIntegralReport_0 = 
      simulation_0.getReportManager().createReport(VolumeIntegralReport.class);

    //FieldHistory fieldHistory_0 = 
     // simulation_0.getMonitorManager().createFieldHistory(userFieldFunction_7);


    userFieldFunction_7.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SCALAR);

    volumeIntegralReport_0.setFieldFunction(userFieldFunction_7);

    volumeIntegralReport_0.setPresentationName("viscous_dissipation");

    volumeIntegralReport_0.printReport();

    volumeIntegralReport_0.printReport();

    volumeIntegralReport_0.getParts().setQuery(null);

    Region region_0 = 
      simulation_0.getRegionManager().getRegion("Region");

    volumeIntegralReport_0.getParts().setObjects(region_0);

    volumeIntegralReport_0.printReport();

// Turbulent Dissipation
	
    UserFieldFunction userFieldFunction_8 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_8.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SCALAR);

    userFieldFunction_8.setPresentationName("turbulent_dissipation");

    userFieldFunction_8.setFunctionName("turbulent_dissipation");

    userFieldFunction_8.setDimensions(Dimensions.Builder().mass(1).length(-1).time(3).build());

    userFieldFunction_8.setDimensions(Dimensions.Builder().mass(1).length(-1).time(-3).build());

    userFieldFunction_8.setDefinition("${Density}*${KwDissipation}");
    
    VolumeIntegralReport volumeIntegralReport_1 = 
      simulation_0.getReportManager().createReport(VolumeIntegralReport.class);

    volumeIntegralReport_1.setFieldFunction(userFieldFunction_8);

    volumeIntegralReport_1.setPresentationName("turbulent_dissipation");

    volumeIntegralReport_1.getParts().setQuery(null);

    volumeIntegralReport_1.getParts().setObjects(region_0);

    volumeIntegralReport_1.printReport();

// Turbulent Viscous Dissipation

    UserFieldFunction userFieldFunction_9 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_9.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SYMMETRIC_TENSOR);

    userFieldFunction_9.setPresentationName("turbulent_stress_tensor");

    userFieldFunction_9.setFunctionName("turbulent_stress_tensor");

    userFieldFunction_9.setDefinition("($$${velocity_gradient_tensor} + $$${velocity_gradient_transpose}) * ${TurbulentViscosity}");

    userFieldFunction_9.setDimensions(Dimensions.Builder().mass(1).length(-1).time(-2).build());

    
    UserFieldFunction userFieldFunction_10 = 
      simulation_0.getFieldFunctionManager().createFieldFunction();

    userFieldFunction_10.getTypeOption().setSelected(FieldFunctionTypeOption.Type.SCALAR);

    userFieldFunction_10.setPresentationName("turbulent_viscous_dissipation");

    userFieldFunction_10.setFunctionName("turbulent_viscous_dissipation");

    userFieldFunction_10.setDimensions(Dimensions.Builder().mass(1).length(-1).time(-3).build());

    userFieldFunction_10.setDefinition("$$${turbulent_stress_tensor}[0,0]*$${grad_vel_u}[0] + $$${turbulent_stress_tensor}[0,1]*$${grad_vel_u}[1] + $$${turbulent_stress_tensor}[0,2]*$${grad_vel_u}[2] + $$${turbulent_stress_tensor}[1,0]*$${grad_vel_v}[0] + $$${turbulent_stress_tensor}[1,1]*$${grad_vel_v}[1] + $$${turbulent_stress_tensor}[1,2]*$${grad_vel_v}[2] + $$${turbulent_stress_tensor}[2,0]*$${grad_vel_w}[0] + $$${turbulent_stress_tensor}[2,1]*$${grad_vel_w}[1] + $$${turbulent_stress_tensor}[2,2]*$${grad_vel_w}[2] ");

    VolumeIntegralReport volumeIntegralReport_2 = 
      simulation_0.getReportManager().createReport(VolumeIntegralReport.class);

    volumeIntegralReport_2.setFieldFunction(userFieldFunction_10);

    volumeIntegralReport_2.setPresentationName("turbulent_viscous_dissipation");

    volumeIntegralReport_2.getParts().setQuery(null);

    volumeIntegralReport_2.getParts().setObjects(region_0);

    volumeIntegralReport_2.printReport();


  }
}
