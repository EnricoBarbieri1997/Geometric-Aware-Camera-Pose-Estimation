{
  # Nix development environment for a Julia project
  description = "Dev shell for Julia project: cylinder-based-homotopy-camera-calibration";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    nix-ld.url = "github:Mic92/nix-ld";
    nix-ld.inputs.nixpkgs.follows = "nixpkgs";
  };

  outputs = { self, nixpkgs, flake-utils, nix-ld }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
      in {
        devShells.default = pkgs.mkShell {
          name = "julia-dev-shell";
          buildInputs = [
            pkgs.julia_111
            pkgs.git
            pkgs.gnumake
            pkgs.python3
            pkgs.gcc
          ];
          shellHook = ''
            echo "Welcome to the Julia dev shell!"
            export JULIA_PROJECT=$PWD
          '';
        };
      }
    );
    nixosConfigurations.nixos = nixpkgs.lib.nixosSystem {
      system = "x86_64-linux";
      modules = [
        nix-ld.nixosModules.nix-ld
        { programs.nix-ld.dev.enable = true; }
      ];
    };
}
