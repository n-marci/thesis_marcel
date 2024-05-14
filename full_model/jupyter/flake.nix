{
  description = "micromamba development flake";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      pkgs = nixpkgs.legacyPackages.x86_64-linux;
      path = "/home/marci/dev/jptest3";
    in {
      devShell.x86_64-linux = (pkgs.buildFHSUserEnv {
        name = "conda";
        targetPkgs = pkgs: (
          with pkgs; [
            micromamba
            jupyter-all
            # gcc
          ]
        );
        profile = ''
          eval "$(micromamba shell hook -s bash)"
          export MAMBA_ROOT_PREFIX=${path}/.mamba

          micromamba create -n calcs jupyter --file requirements.txt -c conda-forge
          micromamba install --yes -f requirements.txt

          micromamba activate calcs
          # micromamba update --yes -f environment.yml
          # jupyter notebook
        '';
      }).env;
    };
}
