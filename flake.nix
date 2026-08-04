{
  description = "rabe pinned build/test toolchain for reproducible golden records";

  # nixos-22.11 pins gfortran 11 + glibc 2.35, the toolchain the committed
  # golden record was produced on; newer compilers/libm drift past rtol=1e-10.
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-22.11";

  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" ];
      forAll = f: nixpkgs.lib.genAttrs systems (s: f nixpkgs.legacyPackages.${s});
      toolchain = pkgs: [
        pkgs.gfortran
        pkgs.cmake
        pkgs.gnumake
        pkgs.pkg-config
        pkgs.git
      ];
      # Python for the golden-record compare; only needed where the test runs.
      pyenv = pkgs: pkgs.python3.withPackages (ps: [ ps.xarray ps.numpy ps.netcdf4 ]);
      mkShell = pkgs: packages:
        pkgs.mkShell {
          packages = packages;
          FC = "${pkgs.gfortran}/bin/gfortran";
        };
    in {
      devShells = forAll (pkgs: {
        default = mkShell pkgs (toolchain pkgs ++ [ (pyenv pkgs) ]);
      });
    };
}
