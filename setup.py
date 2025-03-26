from setuptools import setup

setup(
    name='dsrnascan',
    version='0.1.0',
    py_modules=['dsRNAscan'],
    install_requires=[
        'biopython',
        'numpy',
        'pandas',
        'ViennaRNA'
        # Add additional dependencies as needed (e.g., ViennaRNA bindings)
    ],
    entry_points={
        'console_scripts': [
            'dsrnascan = dsRNAscan:main',
        ],
    },
)
