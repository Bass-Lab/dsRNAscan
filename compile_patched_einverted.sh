#!/bin/bash
# Compile einverted with G-U wobble patch

set -e

echo "Compiling einverted with G-U wobble patch..."

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

# Compile einverted
echo "Compiling einverted..."
cd emboss
make einverted

# Copy to tools directory
echo "Installing patched einverted..."
cp einverted ../../dsrnascan/tools/einverted
chmod +x ../../dsrnascan/tools/einverted

# Also save platform-specific version
PLATFORM=$(uname -s | tr '[:upper:]' '[:lower:]')
ARCH=$(uname -m)

if [[ "$PLATFORM" == "darwin" ]]; then
    if [[ "$ARCH" == "arm64" ]] || [[ "$ARCH" == "aarch64" ]]; then
        cp einverted ../../dsrnascan/tools/einverted_darwin_arm64
    else
        cp einverted ../../dsrnascan/tools/einverted_darwin_x86_64
    fi
elif [[ "$PLATFORM" == "linux" ]]; then
    if [[ "$ARCH" == "aarch64" ]]; then
        cp einverted ../../dsrnascan/tools/einverted_linux_aarch64
    else
        cp einverted ../../dsrnascan/tools/einverted_linux_x86_64
    fi
fi

echo "✓ Successfully compiled einverted with G-U wobble patch!"
echo ""
echo "Testing G-U pairing..."
echo ">test_gu" > test_gu.fa
echo "GGGGGUUUUU" >> test_gu.fa

if ./einverted -sequence test_gu.fa -gap 12 -threshold 10 -match 3 -mismatch -4 2>/dev/null | grep -q "Score: "; then
    echo "✓ G-U pairing detected successfully!"
else
    echo "⚠ Warning: G-U pairing may not be working correctly"
fi

cd ../..
echo "Done! Patched einverted installed in dsrnascan/tools/"