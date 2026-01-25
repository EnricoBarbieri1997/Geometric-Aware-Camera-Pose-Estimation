{
  # Nix development environment for a Julia project
  description = "Dev shell for Julia project: cylinder-based-homotopy-camera-calibration";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
      in {
        devShells.default = pkgs.mkShell {
          name = "julia-dev-shell";
          buildInputs = [
            pkgs.julia
            pkgs.git
            pkgs.gnumake
            pkgs.python3
            pkgs.which
            pkgs.gcc
            pkgs.pkg-config
            pkgs.openssl
            pkgs.zlib
            pkgs.libGL
            pkgs.libGLU
            pkgs.cairo
            pkgs.ffmpeg
            pkgs.glib
            pkgs.xorg.libX11
            pkgs.xorg.libXext
            pkgs.xorg.libXrender
            pkgs.xorg.libXtst
            pkgs.xorg.libXi
            pkgs.xorg.libSM
            pkgs.xorg.libICE
            pkgs.xorg.libxcb
            pkgs.xorg.libXau
            pkgs.xorg.libXdmcp
            pkgs.xorg.libXfixes
            pkgs.xorg.libXrandr
            pkgs.xorg.libXcursor
            pkgs.xorg.libXinerama
            pkgs.xorg.libXcomposite
            pkgs.xorg.libXdamage
            pkgs.xorg.libXxf86vm
            pkgs.xorg.libdrm
            pkgs.xorg.libpciaccess
            pkgs.xorg.libxshmfence
            pkgs.xorg.libxkbfile
            pkgs.xorg.xcbutil
            pkgs.xorg.xcbutilimage
            pkgs.xorg.xcbutilkeysyms
            pkgs.xorg.xcbutilrenderutil
            pkgs.xorg.xcbutilwm
            pkgs.xorg.xcbutilcursor
            pkgs.xorg.xcbutilerrors
            pkgs.xorg.xcbutilxrm
            pkgs.xorg.xcbproto
            pkgs.xorg.xorgproto
            pkgs.xorg.xtrans
            pkgs.xorg.xkeyboardconfig
            pkgs.xorg.xkbcomp
            pkgs.xorg.xmodmap
            pkgs.xorg.xset
            pkgs.xorg.xsetroot
            pkgs.xorg.xrandr
            pkgs.xorg.xev
            pkgs.xorg.xdpyinfo
            pkgs.xorg.xwininfo
            pkgs.xorg.xprop
            pkgs.xorg.xhost
            pkgs.xorg.xauth
            pkgs.xorg.xinit
            pkgs.xorg.xinput
            pkgs.xorg.xkill
            pkgs.xorg.xlsclients
            pkgs.xorg.xmessage
            pkgs.xorg.xpr
            pkgs.xorg.xrefresh
            pkgs.xorg.xwd
            pkgs.xorg.xwud
            pkgs.xorg.xeyes
            pkgs.xorg.xclock
            pkgs.xorg.xlogo
            pkgs.xorg.xcalc
            pkgs.xorg.xclipboard
            pkgs.xorg.xcmsdb
            pkgs.xorg.xcursorgen
            pkgs.xorg.xdriinfo
            pkgs.xorg.xgamma
            pkgs.xorg.xkbprint
            pkgs.xorg.xkbevd
            pkgs.xorg.xkbutils
            pkgs.xorg.xmodmap
            pkgs.xorg.xrdb
            pkgs.xorg.xsetmode
            pkgs.xorg.xsm
            pkgs.xorg.xstdcmap
            pkgs.xorg.xvidtune
            pkgs.xorg.xvinfo
            pkgs.xorg.xwd
            pkgs.xorg.xwininfo
            pkgs.xorg.xwud
          ];
          shellHook = ''
            echo "Welcome to the Julia dev shell!"
            export JULIA_PROJECT=$PWD
          '';
        };
      }
    );
}
