from setuptools import setup

requirements = [
    'numpy', 'pandas', 'h5py',
]

setup(
    name='scRNA_biomarker',
    version='',
    packages=['src', 'src.analysis', 'src.common', 'src.preprocess'],
    url='',
    license='',
    author='NYGC-singlecell-hackathon-team-biomarker',
    author_email='',
    description='',
    install_requires=requirements,
)
