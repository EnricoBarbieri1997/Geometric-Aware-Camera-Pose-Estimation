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
        libPath = pkgs.lib.makeLibraryPath [
          pkgs.stdenv.cc.cc.lib     # libquadmath.so.0, libstdc++.so, etc.
          pkgs.gfortran.cc.lib      # libgfortran.so.*, often needed by JLLs
          
          # Wayland + input
          pkgs.wayland
          pkgs.libxkbcommon

          # EGL/OpenGL (GLVND + Mesa)
          pkgs.libglvnd
          pkgs.mesa

          # Often needed by GLFW/Makie stack
          pkgs.glib
          pkgs.fontconfig
          pkgs.freetype
        ];
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
            export LD_LIBRARY_PATH="${libPath}:$LD_LIBRARY_PATH"

            # Force GLFW to use Wayland (avoid GLX entirely)
            export GLFW_USE_WAYLAND=1

            # Make sure EGL is the chosen GL platform
            export EGL_PLATFORM=wayland

            # Helpful sanity output
            echo "WAYLAND_DISPLAY=$WAYLAND_DISPLAY"
            echo "GLFW_USE_WAYLAND=$GLFW_USE_WAYLAND"
            echo "EGL_PLATFORM=$EGL_PLATFORM"

            exec zsh
          '';
        };
      }
    )
    // {
      nixosConfigurations.nixos = nixpkgs.lib.nixosSystem {
        system = "x86_64-linux";
        modules = [
          nix-ld.nixosModules.nix-ld
          { programs.nix-ld.enable = true; }
        ];
      };
    };
}
