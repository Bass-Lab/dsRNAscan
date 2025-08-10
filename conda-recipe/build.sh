#!/bin/bash
set -e

# Debug info
echo "Building dsRNAscan in $PREFIX"
echo "Source directories:"
ls -la

# Build EMBOSS first
cd emboss_src

# Update config files for arm64 compatibility
curl -L -o config.sub "https://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.sub"
curl -L -o config.guess "https://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.guess"
chmod +x config.sub config.guess

# Configure EMBOSS
./configure --prefix="$PREFIX" --without-x --without-pngdriver

# Build and install EMBOSS
make -j${CPU_COUNT}
make install

# Apply the G-U pairing patch to einverted
cd ../dsrnascan_src
if [ -f "einverted.patch" ]; then
    echo "Applying einverted G-U pairing patch..."
    patch -p0 < einverted.patch -d ../emboss_src/emboss/ || true
    
    # Rebuild einverted with the patch
    cd ../emboss_src/emboss
    make einverted
    cp einverted "$PREFIX/bin/"
    cd ../../dsrnascan_src
fi

# Create tools directory and copy einverted
mkdir -p ./tools
cp "$PREFIX/bin/einverted" ./tools/

# Install the Python package
cd ..
$PYTHON -m pip install ./dsrnascan_src --no-deps -vv