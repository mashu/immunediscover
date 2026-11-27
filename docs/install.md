# Installing

Ready-made programs are published for Linux (Intel/AMD and ARM) and for
macOS (Apple silicon and Intel). Nothing needs to be compiled.

## With the installer

Paste this into a terminal:

```
curl -fsSL https://github.com/mashu/immunediscover/releases/latest/download/install.sh | sh
```

The installer downloads the program for your computer, checks that the
download is intact, puts it into `~/.local/bin` and makes sure that
folder is on your `PATH`. When it has finished, open a new terminal and
run:

```
immunediscover --help
```

### If you already have a program called immunediscover

The installer does not replace or hide it. If another `immunediscover`
(the Julia tool, for example) is found on your `PATH`, this program is
installed as **`immunediscover-rs`**, and the installer tells you so.
Use that name wherever this guide writes `immunediscover`:

```
immunediscover-rs discover blast demux.tsv.gz IGHV.fasta novel.tsv -g V
```

To tell the two apart at any time: this program reports its version as
`immunediscover <version> (Rust)`.

### Options

Options go after `sh -s --`:

```
curl -fsSL https://github.com/mashu/immunediscover/releases/latest/download/install.sh | sh -s -- --bin-dir /opt/tools
```

| Option | What it does |
|---|---|
| `--bin-dir DIR` | Install into this folder instead of `~/.local/bin` (`/usr/local/bin` when run as root). |
| `--name NAME` | Install the program under a name of your choice. |
| `--force` | Install as `immunediscover` even if another program of that name exists. Which of the two then runs depends on the order of your `PATH`. |
| `--version TAG` | Install a particular release, for example `v0.5.0`, instead of the latest. |
| `--no-modify-path` | Do not add the folder to the `PATH` in your shell profile; the installer prints the line to add instead. |
| `--uninstall` | Remove the installed program. |

## From a downloaded file

Each [release](https://github.com/mashu/immunediscover/releases) has one
file per kind of computer:

| File | For |
|---|---|
| `immunediscover-linux-x86_64.tar.gz` | Linux on Intel or AMD processors (most servers and PCs) |
| `immunediscover-linux-aarch64.tar.gz` | Linux on ARM processors |
| `immunediscover-macos-universal.tar.gz` | macOS 11 or newer, Apple silicon and Intel |

Unpack it and run the installer inside, which installs the program that
came with it:

```
tar -xzf immunediscover-linux-x86_64.tar.gz
cd immunediscover-linux-x86_64
./install.sh
```

The program is a single file. You can also copy `immunediscover` from the
unpacked folder to wherever you keep your tools, without the installer.

The Linux programs do not depend on the libraries of the system, so they
run on old and new distributions alike. `SHA256SUMS` in each release
lists the checksums of the files.

## What else is needed

Two commands use other programs, which the installer does not install. It
does check for them, on Linux and macOS alike, and tells you at the end
what it found and what is missing.

| For | You need | One way to get it (Linux and macOS) |
|---|---|---|
| `discover blast` (unless run with `--aligner builtin`) | BLAST+ 2.15 to 2.17: `blastn` and `makeblastdb` on the `PATH` | `conda install -c conda-forge -c bioconda blast=2.16` |
| `search bwa` | `bwa`, and a `bwa index` of the genome | `conda install -c conda-forge -c bioconda bwa` |

All other commands need nothing else. `discover blast --aligner builtin`
and `align` compare reads and alleles with the
[aligner that is part of immunediscover](align.md).

BLAST+ must be one of the versions named: `discover blast` refuses to run
with an older or newer one, because its results depend on the version.
A BLAST+ from another source (Homebrew, the NCBI download, a Linux
package) works as long as its version is in that range.

## Updating and removing

Run the installer again to update to the latest release. To remove the
program:

```
curl -fsSL https://github.com/mashu/immunediscover/releases/latest/download/install.sh | sh -s -- --uninstall
```

## Building from source

With a Rust toolchain (1.82 or newer):

```
cargo build --release
```

The program is then `target/release/immunediscover`.
