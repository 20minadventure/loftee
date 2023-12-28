from setuptools import setup, find_packages


setup(
    name="loftee_annot",
    version="0.1",
    description="",
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'annotate_vcf = analysis.cmd.split_vep:annotate_vcf',
            'rare_variants_table = analysis.cmd.split_vep:rare_variants_table',
            'install_vep = preprocessing.install_vep:main',
        ]
    },
)
