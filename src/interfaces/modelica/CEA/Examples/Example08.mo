within CEA.Examples;
model Example08
  extends Modelica.Icons.Example;
  import CEA.Interfaces.*;

  CEA_Problem prob = CEA_Problem();

algorithm
  Modelica.Utilities.Streams.print("hello");
  Modelica.Utilities.Streams.print("size = " + String(ffi_cea_sizeof(prob)));
end Example08;
