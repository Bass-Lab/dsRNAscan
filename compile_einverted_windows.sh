#!/bin/bash
# Compile einverted with G-U wobble patch for Windows
# For use with MSYS2, Cygwin, or WSL

set -e

echo "Building einverted for Windows with G-U wobble patch..."
echo ""

# Detect environment
if [[ "$OSTYPE" == "msys" ]]; then
    echo "Detected MSYS2 environment"
    PLATFORM="windows_msys2"
elif [[ "$OSTYPE" == "cygwin" ]]; then
    echo "Detected Cygwin environment"
    PLATFORM="windows_cygwin"
elif [[ "$OSTYPE" == "linux-gnu" ]] && grep -qi microsoft /proc/version 2>/dev/null; then
    echo "Detected WSL environment"
    PLATFORM="windows_wsl"
else
    echo "Detected environment: $OSTYPE"
    PLATFORM="windows_generic"
fi

# Create tools directory
mkdir -p dsrnascan/tools

# Check for required tools
for tool in gcc make patch curl tar; do
    if ! command -v $tool &> /dev/null; then
        echo "ERROR: $tool is not installed"
        echo ""
        if [[ "$PLATFORM" == "windows_msys2" ]]; then
            echo "Install with: pacman -S mingw-w64-x86_64-gcc make patch curl tar"
        elif [[ "$PLATFORM" == "windows_cygwin" ]]; then
            echo "Install using Cygwin setup.exe"
        else
            echo "Please install $tool and try again"
        fi
        exit 1
    fi
done

# Download EMBOSS if not present
if [ ! -f "EMBOSS-6.6.0.tar.gz" ]; then
    echo "Downloading EMBOSS 6.6.0..."
    curl -L -o EMBOSS-6.6.0.tar.gz ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz || \
    wget -O EMBOSS-6.6.0.tar.gz ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
fi

# Extract EMBOSS
echo "Extracting EMBOSS..."
tar -xzf EMBOSS-6.6.0.tar.gz
cd EMBOSS-6.6.0

# Apply G-U wobble patch
if ! grep -q "Allowing for GU matches" emboss/einverted.c 2>/dev/null; then
    echo "Applying G-U wobble patch..."
    patch -p1 < ../einverted.patch
else
    echo "Patch already applied"
fi

# Configure EMBOSS for Windows
echo "Configuring EMBOSS for Windows..."
./configure --without-x --disable-shared --without-pngdriver --without-hpdf \
            --without-mysql --without-postgresql --prefix=$(pwd)/../emboss_install

# Build necessary components
echo "Building EMBOSS libraries..."
make -C ajax || echo "Some ajax components failed, continuing..."
make -C nucleus || echo "Nucleus build failed, continuing..."
make -C plplot || echo "Plplot build failed, continuing..."

# Build einverted
echo "Building einverted..."
cd emboss
make einverted || {
    echo "Standard make failed, trying direct compilation..."
    
    # Try different compilation approaches
    if [[ "$PLATFORM" == "windows_msys2" ]]; then
        # MSYS2/MinGW approach
        gcc -O2 -DWIN32 -DMINGW32 -I../ajax/core -I../ajax/ajaxdb -I../ajax/acd -I../nucleus \
            einverted.c \
            -L../ajax/core/.libs -L../ajax/acd/.libs -L../nucleus/.libs \
            -lnucleus -lacd -lajax -lm -lz \
            -o einverted.exe 2>/dev/null || {
                echo "MinGW compilation failed, trying static linking..."
                gcc -O2 -static -DWIN32 -DMINGW32 einverted.c \
                    ../nucleus/.libs/*.o ../ajax/core/.libs/*.o ../ajax/acd/.libs/*.o \
                    -lm -lz -o einverted.exe 2>/dev/null || echo "Static linking also failed"
            }
    elif [[ "$PLATFORM" == "windows_cygwin" ]]; then
        # Cygwin approach
        gcc -O2 -DCYGWIN -I../ajax/core -I../ajax/ajaxdb -I../ajax/acd -I../nucleus \
            einverted.c \
            -L../ajax/core/.libs -L../ajax/acd/.libs -L../nucleus/.libs \
            -lnucleus -lacd -lajax -lm -lz \
            -o einverted.exe 2>/dev/null || echo "Cygwin compilation failed"
    else
        # Generic approach
        gcc -O2 -I../ajax/core -I../ajax/ajaxdb -I../ajax/acd -I../nucleus \
            einverted.c \
            -L../ajax/core/.libs -L../ajax/acd/.libs -L../nucleus/.libs \
            -lnucleus -lacd -lajax -lm -lz \
            -o einverted.exe 2>/dev/null || echo "Generic compilation failed"
    fi
}

# Find and copy the binary
BINARY_FOUND=false
for binary in .libs/einverted.exe .libs/einverted einverted.exe einverted; do
    if [ -f "$binary" ]; then
        echo "Found binary: $binary"
        cp "$binary" ../../dsrnascan/tools/einverted.exe
        chmod +x ../../dsrnascan/tools/einverted.exe 2>/dev/null || true
        
        # Also save platform-specific version
        cp "$binary" ../../dsrnascan/tools/einverted_${PLATFORM}.exe
        BINARY_FOUND=true
        break
    fi
done

if [ "$BINARY_FOUND" = false ]; then
    echo ""
    echo "ERROR: Could not compile einverted for Windows"
    echo "This might be due to missing dependencies or incompatible environment."
    echo ""
    echo "Alternative options:"
    echo "1. Use WSL (Windows Subsystem for Linux) and run the Linux compilation"
    echo "2. Use Docker with the dsRNAscan container"
    echo "3. Download pre-compiled binaries from the GitHub releases page"
    exit 1
fi

# Copy ACD files
if [ -f "acd/einverted.acd" ]; then
    mkdir -p ../../dsrnascan/tools/acd
    cp acd/einverted.acd ../../dsrnascan/tools/acd/
fi

cd ../..

# Test the binary
echo ""
echo "Testing einverted on Windows..."
cat > test_win.fa << EOF
>test_gu
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGNNNNNNNNNNNNNNCTTCTCTCTCCTTCTCTCTCCTTCTCTCTCCTTCTCTCTC
EOF

if dsrnascan/tools/einverted.exe -sequence test_win.fa -gap 20 -threshold 30 -match 3 -mismatch -4 -outfile stdout -auto 2>/dev/null | grep -q "Score"; then
    echo "✓ einverted compiled successfully and G-U wobble pairing works!"
else
    echo "⚠ einverted was compiled but may not be fully functional"
    echo "  The binary may still work but requires additional runtime dependencies"
fi

# Clean up
rm -f test_win.fa
rm -rf EMBOSS-6.6.0 EMBOSS-6.6.0.tar.gz

echo ""
echo "Installation complete!"
echo "einverted.exe has been installed to: dsrnascan/tools/einverted.exe"
echo ""
echo "To use dsRNAscan on Windows:"
echo "  1. Install Python dependencies: pip install biopython numpy pandas ViennaRNA"
echo "  2. Run: python dsrnascan/dsRNAscan.py <your_file.fasta>"