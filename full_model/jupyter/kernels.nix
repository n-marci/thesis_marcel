{pkgs, ...}: {
  kernel.python.minimal = {
    enable = true;
    # extraPackages = ps: [ ps.dulwich ];
    # requiredRuntimePackages = [ pkgs.poetry ];
  };
  kernel.python.python_plus_modules = {
    enable = true;
    projectDir = ./python_plus_modules;
    requiredRuntimePackages = [ pkgs.poetry ];
  };
}
