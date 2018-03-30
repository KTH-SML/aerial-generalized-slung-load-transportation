### Slung load transportation 

The notebooks, described next, are meant to be read in conjunction with the companion pdf document "aerial generalized slung load transportation.pdf".

The model for the system, composed of two rigid bodies (UAV and rod-like manipulator), is described in the notebook "vector_field.nb"  (Section IV).

Two cases are considered:
- torque input on manipulator is available -- UAV-manipulator system
- torque input on manipulator is not available -- UAV-slung-manipulator system

In the notebook
- "uav_manipulator_transformation.nb", we present the transformation for the UAV-manipulator system that decomposes the system into two decoupled systems: one thrust-propelled system, and one unit vector double integrator system  (Section VI) 
- "uav_slung_manipulator_transformation.nb", we present the transformation for the UAV-slung-manipulator system that decomposes the system into one thrust-propelled system, cascaded after one unit vector double integrator system  (Section VII) 

In the notebook
- "uav_manipulator_control.nb", we present the specific controllers for the transformed UAV-manipulator system (Section VI)
- "uav_slung_manipulator_control.nb", we present the specific controllers for the transformed UAV-slung-manipulator system (Section VII)
