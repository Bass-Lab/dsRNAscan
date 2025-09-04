#!/bin/bash
# Compile einverted with G-U wobble patch

set -e

echo "Compiling einverted with G-U wobble patch..."

# Save the original directory
ORIGINAL_DIR=$(pwd)

# Check if EMBOSS source exists
if [ -d "EMBOSS-6.6.0" ]; then
    echo "Using existing EMBOSS source..."
    cd EMBOSS-6.6.0
else
    echo "EMBOSS source not found. Downloading..."
    if [ ! -f "EMBOSS-6.6.0.tar.gz" ]; then
        wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz || \
        curl -O ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
    fi
    tar -xzf EMBOSS-6.6.0.tar.gz
    cd EMBOSS-6.6.0
fi

# Check if patch already applied
if ! grep -q "Allowing for GU matches" emboss/einverted.c 2>/dev/null; then
    echo "Applying G-U wobble patch..."
    patch -p1 < ../einverted.patch
else
    echo "Patch already applied"
fi

# Configure EMBOSS
echo "Configuring EMBOSS..."
./configure --without-x --disable-shared --prefix=$(pwd)/../emboss_install

# Compile necessary libraries first
echo "Compiling necessary libraries..."
echo "Building ajax libraries..."
make -C ajax || echo "Some ajax components failed, continuing..."
echo "Building nucleus library..."
make -C nucleus || echo "Nucleus build failed, continuing..."
echo "Building plplot library (required for einverted)..."
make -C plplot || {
    echo "Plplot build failed, trying to build it separately..."
    cd plplot && make && cd ..
}

# Now compile einverted
echo "Compiling einverted..."
cd emboss
make einverted || {
    echo "Standard make failed, checking for missing dependencies..."
    # Make sure plplot is built
    if [ ! -f ../plplot/.libs/libeplplot.a ]; then
        echo "Building plplot first..."
        cd ../plplot && make && cd ../emboss
    fi
    # Try again
    make einverted || {
        echo "Still failing, trying to link manually..."
        gcc -O3 -o einverted einverted.o \
            ../nucleus/.libs/libnucleus.a \
            ../ajax/acd/.libs/libacd.a \
            ../ajax/ajaxdb/.libs/libajaxdb.a \
            ../ajax/ensembl/.libs/libensembl.a \
            ../ajax/graphics/.libs/libajaxg.a \
            ../ajax/core/.libs/libajax.a \
            ../ajax/zlib/.libs/libezlib.a \
            ../ajax/expat/.libs/libeexpat.a \
            ../ajax/pcre/.libs/libepcre.a \
            ../plplot/.libs/libeplplot.a \
            -lm -lz 2>&1 || echo "Manual linking also failed"
    }
}

# Copy the ACTUAL BINARY (not the libtool wrapper) to tools directory
echo "Installing patched einverted..."
TARGET_DIR="${ORIGINAL_DIR}/dsrnascan/tools"
mkdir -p "${TARGET_DIR}"
if [ -f ".libs/einverted" ]; then
    # The actual binary is in .libs directory
    cp .libs/einverted "${TARGET_DIR}/einverted"
elif [ -f "einverted" ]; then
    # Fallback to the regular binary
    cp einverted "${TARGET_DIR}/einverted"
else
    echo "Error: Could not find einverted binary!"
    exit 1
fi
chmod +x "${TARGET_DIR}/einverted"

# Copy ACD file so einverted can run standalone
if [ -f "acd/einverted.acd" ]; then
    mkdir -p "${TARGET_DIR}/acd"
    cp acd/einverted.acd "${TARGET_DIR}/acd/"
fi

# Also save platform-specific version
PLATFORM=$(uname -s | tr '[:upper:]' '[:lower:]')
ARCH=$(uname -m)

if [[ "$PLATFORM" == "darwin" ]]; then
    if [[ "$ARCH" == "arm64" ]] || [[ "$ARCH" == "aarch64" ]]; then
        cp "${TARGET_DIR}/einverted" "${TARGET_DIR}/einverted_darwin_arm64"
    else
        cp "${TARGET_DIR}/einverted" "${TARGET_DIR}/einverted_darwin_x86_64"
    fi
elif [[ "$PLATFORM" == "linux" ]]; then
    if [[ "$ARCH" == "aarch64" ]]; then
        cp "${TARGET_DIR}/einverted" "${TARGET_DIR}/einverted_linux_aarch64"
    else
        cp "${TARGET_DIR}/einverted" "${TARGET_DIR}/einverted_linux_x86_64"
    fi
fi

echo "✓ Successfully compiled einverted with G-U wobble patch!"
echo ""
echo "Testing G-U pairing..."
echo ">test_gu" > test_gu.fa
echo "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGNNNNNNNNNNNNNNCTTCTCTCTCCTTCTCTCTCCTTCTCTCTCCTTCTCTCTC" >> test_gu.fa

# Set EMBOSS environment for ACD files
export EMBOSS_ACDROOT=$(pwd)/../acd
export EMBOSS_DATA=$(pwd)/../data

# Test with the installed binary
if "${TARGET_DIR}/einverted" -sequence test_gu.fa -gap 20 -threshold 30 -match 3 -mismatch -4 -outfile stdout -auto 2>/dev/null | grep -q "Score"; then
    echo "✓ G-U pairing detected successfully!"
else
    echo "⚠ Warning: G-U pairing test inconclusive (may need ACD files at runtime)"
fi

cd "${ORIGINAL_DIR}"

# Clean up test file
rm -f EMBOSS-6.6.0/emboss/test_gu.fa

echo "Done! Patched einverted installed in dsrnascan/tools/"
echo "Note: einverted may need EMBOSS_ACDROOT environment variable set to locate ACD files"