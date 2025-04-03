#!/bin/bash
set -e

# Set PREFIX to a directory within the dsRNAscan folder if not already defined
if [ -z "$PREFIX" ]; then
    export PREFIX="$(pwd)/emboss_install"
fi

# Determine number of physical CPUs (fallback to sysctl if nproc not available)
if command -v nproc >/dev/null 2>&1; then
    CORES=$(nproc)
else
    CORES=$(sysctl -n hw.physicalcpu)
fi

# Define the EMBOSS tarball URL and directory name
EMBOSS_URL="ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz"
EMBOSS_DIR="EMBOSS-6.6.0"

# Remove any existing EMBOSS directory to ensure a fresh build
if [ -d "$EMBOSS_DIR" ]; then
    rm -rf "$EMBOSS_DIR"
fi

# Download the EMBOSS tarball if it doesn't exist
if [ ! -f "EMBOSS-6.6.0.tar.gz" ]; then
    echo "Downloading EMBOSS-6.6.0..."
    curl -O "$EMBOSS_URL"
fi

# Extract the tarball
tar -xzf EMBOSS-6.6.0.tar.gz

# Change into the EMBOSS source directory
cd "$EMBOSS_DIR"

# Clean previous builds if any (if config.cache exists, remove it and do make clean)
if [ -f config.cache ]; then
    rm config.cache
    make clean
fi

# Update config.sub and config.guess for arm64 compatibility
curl -L -o config.sub "https://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.sub"
curl -L -o config.guess "https://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.guess"
chmod +x config.sub config.guess

# Configure EMBOSS with the installation prefix.
# Optionally Disable X11 and plplot to avoid unwanted dependencies using --without-x --without-plplot
./configure --prefix="$PREFIX" 

# Build and install using the determined number of cores
make -j"$CORES"
make install SHAREDIR="$PREFIX/share"

# Return to the repository root directory
cd ..

# Copy the einverted binary to a location within the dsRNAscan package
mkdir -p ./tools
# Assuming that the einverted binary was installed into $PREFIX/bin, copy it:
cp "$PREFIX/bin/einverted" ./tools/

# Install your Python package using pip
$PYTHON -m pip install . --no-deps -vv