Python Environment
------------------
This protein conformational clustering workflow requires a scientific python
installation. We use python 2.7.3. Some additional packages are required.
On a debian based system you can get all of the required packages with

```
sudo apt-get install python-dev python-numpy python-scipy python-setuptools cython
```

If you don't have sudo access or a package manager, we've had good experiences
with the Enthought python distribution (https://www.enthought.com/products/canopy/)
which contains all of these packages preinstalled.

To read the trajectory files from their binary format on disk, you'll need to install
a more specialized package, mdtraj.

```
git clone git://github.com/rmcgibbo/mdtraj.git
cd mdtraj
sudo python setup.py install
```

Getting the data
----------------
Download the file xtc.tar.bz (2.9 gb) from this site: http://purl.stanford.edu/bd829sf1034

```
wget https://stacks.stanford.edu/file/druid:bd829sf1034/xtc.tar.bz
tar -xjvf xtc.tar.bz
```

(Download it directly into this project's directory, since the paths are hardcoded in the cluster.py
script)

Running the code
----------------
This repository contains pure-python implementations of some of the algorithms
in MSMBuilder. We figured the pure python code would be easier to hack than
optimized C. To get started, try running the file cluster.py

```
python cluster.py
```

And then open it up and take a look. It's opening up some protein conformations
and computing all of the pairwise distances between the conformations, using the
RMSD distance metric (implemented in rmsd.py), and then using those distances
as input for a hierarchical clustering.
