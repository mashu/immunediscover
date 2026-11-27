#!/bin/sh
# Installer of immunediscover for Linux and macOS.
#
#   From a downloaded release:   ./install.sh
#   Straight from the web:       curl -fsSL https://github.com/mashu/immunediscover/releases/latest/download/install.sh | sh
#
# Puts the program into a folder on your PATH. If another program called
# "immunediscover" is already installed (the Julia tool, for example), this
# one is installed as "immunediscover-rs" and the other is left alone.

set -eu

REPO="mashu/immunediscover"
PROGRAM="immunediscover"
FALLBACK_NAME="immunediscover-rs"

name=""
bin_dir=""
version="latest"
force=0
modify_path=1
uninstall=0

usage() {
    cat <<USAGE
Usage: install.sh [options]

  --bin-dir DIR      Folder to install into (default: ~/.local/bin,
                     or /usr/local/bin when run as root)
  --name NAME        Install the program under this name
  --force            Install as "$PROGRAM" even if another program of that
                     name is already on the PATH
  --version TAG      Release to download, for example v0.5.0 (default: latest)
  --no-modify-path   Do not add the folder to the PATH in your shell profile
  --uninstall        Remove the installed program
  -h, --help         Show this help
USAGE
}

say() { printf '%s\n' "$*"; }
fail() { printf 'install.sh: %s\n' "$*" >&2; exit 1; }
have() { command -v "$1" >/dev/null 2>&1; }

while [ $# -gt 0 ]; do
    case "$1" in
        --bin-dir) [ $# -ge 2 ] || fail "--bin-dir needs a folder"; bin_dir="$2"; shift 2 ;;
        --name) [ $# -ge 2 ] || fail "--name needs a name"; name="$2"; shift 2 ;;
        --version) [ $# -ge 2 ] || fail "--version needs a tag"; version="$2"; shift 2 ;;
        --force) force=1; shift ;;
        --no-modify-path) modify_path=0; shift ;;
        --uninstall) uninstall=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) usage >&2; fail "unknown option: $1" ;;
    esac
done

if [ -z "$bin_dir" ]; then
    if [ "$(id -u)" = "0" ]; then
        bin_dir="/usr/local/bin"
    else
        [ -n "${HOME:-}" ] || fail "HOME is not set; give --bin-dir"
        bin_dir="$HOME/.local/bin"
    fi
fi

# Is the program at this path ours? Ours says so in its version.
is_ours() {
    "$1" --version 2>/dev/null | grep -q "(Rust)"
}

if [ "$uninstall" = "1" ]; then
    removed=0
    for candidate in "${name:-$PROGRAM}" "$FALLBACK_NAME"; do
        path="$bin_dir/$candidate"
        if [ -f "$path" ] && is_ours "$path"; then
            rm -f "$path"
            say "Removed $path"
            removed=1
        fi
    done
    [ "$removed" = "1" ] || say "Nothing to remove in $bin_dir"
    exit 0
fi

# Which release file fits this computer?
platform() {
    os="$(uname -s)"
    machine="$(uname -m)"
    case "$os" in
        Linux)
            case "$machine" in
                x86_64|amd64) say "linux-x86_64" ;;
                aarch64|arm64) say "linux-aarch64" ;;
                *) fail "no release for Linux on $machine; build from source with cargo" ;;
            esac ;;
        Darwin) say "macos-universal" ;;
        *) fail "no release for $os; build from source with cargo" ;;
    esac
}

download() {
    # download URL FILE
    if have curl; then
        curl -fsSL "$1" -o "$2"
    elif have wget; then
        wget -q "$1" -O "$2"
    else
        fail "curl or wget is needed to download the release"
    fi
}

sha256_of() {
    if have sha256sum; then
        sha256sum "$1" | cut -d' ' -f1
    elif have shasum; then
        shasum -a 256 "$1" | cut -d' ' -f1
    else
        say ""
    fi
}

work=""
cleanup() { [ -z "$work" ] || rm -rf "$work"; }
trap cleanup EXIT

# The program to install: the one next to this script (a downloaded and
# unpacked release), or else the release for this computer from the web.
source_binary=""
case "$0" in
    */install.sh|install.sh)
        script_dir="$(cd "$(dirname "$0")" && pwd)"
        if [ -f "$script_dir/$PROGRAM" ]; then
            source_binary="$script_dir/$PROGRAM"
        fi ;;
esac

if [ -z "$source_binary" ]; then
    asset="$PROGRAM-$(platform).tar.gz"
    if [ "$version" = "latest" ]; then
        base="https://github.com/$REPO/releases/latest/download"
    else
        base="https://github.com/$REPO/releases/download/$version"
    fi
    work="$(mktemp -d)"
    say "Downloading $asset ($version)"
    download "$base/$asset" "$work/$asset" || fail "could not download $base/$asset"
    if download "$base/SHA256SUMS" "$work/SHA256SUMS" 2>/dev/null; then
        expected="$(grep " $asset\$" "$work/SHA256SUMS" | cut -d' ' -f1 || true)"
        actual="$(sha256_of "$work/$asset")"
        if [ -n "$expected" ] && [ -n "$actual" ]; then
            [ "$expected" = "$actual" ] || fail "the download is damaged (checksum mismatch); try again"
        fi
    fi
    tar -xzf "$work/$asset" -C "$work"
    source_binary="$(find "$work" -type f -name "$PROGRAM" | head -n 1)"
    [ -n "$source_binary" ] || fail "the release file does not contain $PROGRAM"
fi

[ -x "$source_binary" ] || chmod +x "$source_binary"
# A file downloaded with a browser on macOS is quarantined; lift that.
if [ "$(uname -s)" = "Darwin" ] && have xattr; then
    xattr -d com.apple.quarantine "$source_binary" 2>/dev/null || true
fi
"$source_binary" --version >/dev/null 2>&1 || fail "the program does not run on this computer"

# Under which name? Another "immunediscover" on the PATH is not replaced
# and not shadowed.
if [ -z "$name" ]; then
    name="$PROGRAM"
    existing="$(command -v "$PROGRAM" 2>/dev/null || true)"
    if [ -n "$existing" ] && [ "$existing" != "$bin_dir/$PROGRAM" ] && ! is_ours "$existing"; then
        if [ "$force" = "1" ]; then
            say "Note: another $PROGRAM is installed at $existing."
            say "      Which of the two runs depends on the order of your PATH."
        else
            name="$FALLBACK_NAME"
            say "Another program called $PROGRAM is already installed:"
            say "    $existing"
            say "It is left as it is. This one is installed as $name."
            say "(Use --force to install it as $PROGRAM anyway, or --name to choose a name.)"
        fi
    fi
fi

mkdir -p "$bin_dir" 2>/dev/null || fail "cannot create $bin_dir; give --bin-dir or run with sudo"
[ -w "$bin_dir" ] || fail "cannot write to $bin_dir; give --bin-dir or run with sudo"
# Copy next to the destination and move into place, so that a running
# program is never half overwritten.
cp "$source_binary" "$bin_dir/.$name.new"
chmod 755 "$bin_dir/.$name.new"
mv -f "$bin_dir/.$name.new" "$bin_dir/$name"

installed="$("$bin_dir/$name" --version)"
say "Installed $installed as $bin_dir/$name"

# Is the folder on the PATH?
case ":$PATH:" in
    *":$bin_dir:"*) on_path=1 ;;
    *) on_path=0 ;;
esac
if [ "$on_path" = "0" ]; then
    line="export PATH=\"$bin_dir:\$PATH\""
    profile=""
    if [ "$modify_path" = "1" ] && [ -n "${HOME:-}" ]; then
        case "$(basename "${SHELL:-sh}")" in
            zsh) profile="$HOME/.zshrc" ;;
            bash)
                if [ "$(uname -s)" = "Darwin" ]; then
                    profile="$HOME/.bash_profile"
                else
                    profile="$HOME/.bashrc"
                fi ;;
            *) profile="$HOME/.profile" ;;
        esac
    fi
    if [ -n "$profile" ]; then
        if ! grep -qsF "$line" "$profile"; then
            printf '\n# Added by the immunediscover installer\n%s\n' "$line" >> "$profile"
        fi
        say "$bin_dir was added to your PATH in $profile."
        say "Open a new terminal, then run: $name --help"
    else
        say "$bin_dir is not on your PATH. Add this line to your shell profile:"
        say "    $line"
    fi
else
    say "Run: $name --help"
fi

if ! have blastn; then
    say ""
    say "Note: 'discover blast' needs BLAST+ 2.15 to 2.17 (blastn, makeblastdb), which was"
    say "not found. One way to get it: conda install -c bioconda blast=2.16"
fi
