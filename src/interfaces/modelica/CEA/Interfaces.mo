within CEA;
package Interfaces
  extends Modelica.Icons.Package;

  class CEA_Problem
    extends ExternalObject;

    function constructor
      output CEA_Problem prob;
    external "C" prob = _ffi_cea_new_problem()
        annotation(Include = "#include \"ffi_cea_modelica.h\"",
                   Library = {"cea_modelica", "-lc++", "-lc++abi", "-lgfortran"},
                   IncludeDirectory = "modelica://CEA/Resources/Include",
                   LibraryDirectory = "modelica://CEA/Resources/Library");
    end constructor;

    function destructor
      input  CEA_Problem prob;
    external "C" _ffi_cea_del_problem(prob);
    end destructor;

  end CEA_Problem;

  function ffi_cea_sizeof
    input  CEA_Problem prob;
    output Integer n;
  external "C" n = _ffi_cea_sizeof_int(prob);
  end ffi_cea_sizeof;

end Interfaces;
