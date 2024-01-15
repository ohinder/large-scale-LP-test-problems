let
  unstable = import
    (fetchTarball "https://nixos.org/channels/nixos-unstable/nixexprs.tar.xz")
    { };
in { pkgs ? import <nixpkgs> { } }:
  with pkgs;
  mkShell { buildInputs = [ unstable.julia marksman ]; }
