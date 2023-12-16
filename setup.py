from setuptools import setup, find_packages


setup(
    name="loftee_annot",
    version="0.1",
    description="",
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'split_vep = analysis.split_vep:main',
            'install_vep = preprocessing.install_vep:main',
        ]
    },
)
