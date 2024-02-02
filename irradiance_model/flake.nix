{
  description = "micromamba development flake";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      pkgs = nixpkgs.legacyPackages.x86_64-linux;
      path = "/home/marci/dev/ma/irradiance_model";
    in {
      devShell.x86_64-linux = (pkgs.buildFHSUserEnv {
        name = "conda";
        targetPkgs = pkgs: (
          with pkgs; [
            micromamba
            # gcc
            taplo
            yaml-language-server
          ]
        );
        profile = ''
          eval "$(micromamba shell hook -s bash)"
          export MAMBA_ROOT_PREFIX=${path}/.mamba

          # micromamba create -n env
          # micromamba install --yes -f environment.yml

          micromamba activate env
          # micromamba update --yes -f environment.yml
        '';
      }).env;
    };
}
