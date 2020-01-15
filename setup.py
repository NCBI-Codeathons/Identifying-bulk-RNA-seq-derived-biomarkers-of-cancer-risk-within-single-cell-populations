from setuptools import setup

requirements = [
    'numpy', 'pandas', 'h5py', 'dash'
]

setup(
    name='scRNA-seq-bulk-biomarker',
    version='',
    packages=['src', 'src.analysis', 'src.common', 'src.preprocess', 'src.visualization'],
    url='',
    license='',
    author='NYGC-singlecell-hackathon-team-biomarker',
    author_email='',
    description='',
    install_requires=requirements,
)
