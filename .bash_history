cd ..
ll
rm -rf Bec3D
ll
module load git/2.9.5
module list
git clone git@github.com:DenMuj/Bec3D.git
cd Bec3D/
git init
cd MojOutput/
ll
cat 637819.paradox.ipb.ac.rs.out
rm *.err
rm *.out
ll
qsub /home/denis/Bec3D/imag.pbs 
cat imag3d-rms.txt
cd bec3d/
ll
module load git/2.9.5 
module avail
module load gnu/9.3.0 
make bec-gp-rot-3d-th compiler=gnu
ll
cd MojOutput/
qsub /home/denis/bec3d/imag.pbs 
ll
cat 637834.paradox.ipb.ac.rs.err
module avail
module load git/2.9.5
git pull
rm *.err
rm *.out
qsub /home/denis/bec3d/imag.pbs 
ll
cat 637835.paradox.ipb.ac.rs.err
rm *.out
rm *.err
cd ..
ll
rm *.o
ll
rm bec-gp-rot-3d-th
ll
make bec-gp-rot-3d-th 
ll
cd MojOutput/
qsub /home/denis/bec3d/imag.pbs 
ll
cat 637836.paradox.ipb.ac.rs.err
cd bec3d/
module load git/2.9.5 
git pull
rm /input/real3d-input
rm ./input/real3d-input
git pull
cd input
ll
cd ..
cd MojOutput/
rm *.out
rm *.err
qsub /home/denis/bec3d/job.pbs 
ll
cat 637843.paradox.ipb.ac.rs.err
rm *.out
rm *.err
git pull
qsub /home/denis/bec3d/job.pbs 
ll
cat 637844.paradox.ipb.ac.rs.err
module avail
git pull
rm *.out
rm *.err
qsub /home/denis/bec3d/job.pbs 
ll
cat 637845.paradox.ipb.ac.rs.err
pip install scipy
rm *.out
rm *.err
git pull
ll
qsub /home/denis/bec3d/job.pbs 
ll
rm *.err
rm *.out
git pull
qsub /home/denis/bec3d/job.pbs 
ll
showq
ll
cat real3d-rms.txt
showq
cat real3d-rms.txt
qdel 637847
git pull
rm *.out
rm *.err
qsub /home/denis/bec3d/job.pbs 
ll
cat real3d-rms.txt 
ll
cat 637848.paradox.ipb.ac.rs.out
cd ..
git pull
rm ./input/real3d-input 
git pull
cd MojOutput/
rm *.out
rm *.err
qsub /home/denis/bec3d/job.pbs 
ll
cat real3d-rms.txt 
showq
qdel 637849
cd bec3d/
module load git/2.9.5 
module load python/3.6.5 
pip install matplotlib
git pull
cd MojOutput/
rm *.err
rm *.out
ll
git pull
qsub /home/denis/bec3d/job.pbs 
ll
cd bec3d/
qstat -au denis
showq
cat MojOutput/real3d-rms.txt 
nproc
mpstat
showq
top
mpstat
cat MojOutput/real3d-rms.txt 
cd MojOutput/
cd bec3d/
cd MojOutput/
cat real3d-rms.txt
ll
cat 637842.paradox.ipb.ac.rs.out
git pull
git add .
git commit -m "FPS"
git push
cdc bec3d/
cd bec3d/
cd MojOutput/
ll
cat 637850.paradox.ipb.ac.rs.out
module load git/2.9.5 
git pull
showq
qsub /home/denis/bec3d/job.pbs 
ll
637854.paradox.ipb.ac.rs.err
cat 637854.paradox.ipb.ac.rs.err
rm *.err
rm *.out
git pull
qsub /home/denis/bec3d/job.pbs 
ll
showq
ll
cat real3d-rms.txt 
cd bec3d/
cd MojOutput/
showq
qdel 637855
showq
ll
rm *.err
rm *.out
git pull
qsub /home/denis/bec3d/job.pbs 
ll
cat real3d-rms.txt 
cd bec3d/
cd MojOutput/
ll
cat real3d-rms.txt 
showq
cd bec3d/
cd MojOutput/
ll
cat real3d-rms.txt 
showq
cd bec3d/
cd MojOutput/
ll
cat real3d-rms.txt 
showq
cd bec3d/
cd MojOutput/
cat real3d-rms.txt 
showq
ll
module avail
module load python/3.6.5
python --version
module purge
module load anaconda/3-5.1.0
module list
jupyter lab
module load anaconda/3-5.1.0
conda search
module list
module load anaconda/2-5.1.0 
conda search SciPy
module load git/2.9.5 
module avail
git pull
cd bec3d/
git pull
qsub anac.pbs
ll
showq
ll
showq
ll
qdel 637857
rm anac.pbs 
showq
cd bec3d/
cd MojOutput/
ll
rm *.err
rm *.out
ll
module load git/2.9.5 
git add . 
git commit -m "Moje"
car real3d-rms.txt 
cat real3d-rms.txt 
cd bec3d/
cd MojOutput/
cd plots/
ll
fer my_plot.png
feh my_plot.png
jp2a my_plot.png
cd bec3d/MojOutput/
ll
cd bec3d/MojOutput/
ll
cat real3d-rms.txt 
showq
cd bec3d/
cd MojOutput/
ll
cat real3d-rms.txt 
showq
module avail
module load git/2.9.5 
cd ..
git pull
cd MojOutput/
git pull
qsub anac.pbs
qsub /home/denis/bec3d/anac.pbs
ll
cat 637858.paradox.ipb.ac.rs.err
git pull
rm *.err
rm *.out
qsub /home/denis/bec3d/anac.pbs
ll
showq
ll
qdel 637859
showq
git pull
ll
rm *.err
rm *.out
qsub /home/denis/bec3d/anac.pbs
ll
cat real3d-rms.txt 
ll
showq
ll
l
ll
showq
qdel 637860
ll
rm *.err
rm *.out
git pull
qsub /home/denis/bec3d/anac.pbs
ll
top
ll
cd ..
ll
rm *.err
rm *.out
cd MojOutput/
ll
showq
ll
cat 637861.paradox.ipb.ac.rs.err
ll
rm *.err
rm *.out
git pull
rm /home/denis/bec3d/ana.py 
git pull
qsub /home/denis/bec3d/anac.pbs
ll
cat 637862.paradox.ipb.ac.rs.err
git pull
rm /home/denis/bec3d/ana.py 
git pull
rm *.err
rm *.out
qsub /home/denis/bec3d/anac.pbs
ll
cd plots/
ll
cd ..
cat /home/denis/bec3d/input/real3d-input 
cd bec3d/MojOutput/
ll
car real3d-rms.txt 
cat real3d-rms.txt 
cat /home/denis/bec3d/input/real3d-input 
cd bec3d/MojOutput/
ll
cat 637856.paradox.ipb.ac.rs.out
showq
module avail
showq
ll
cd bec3d/
cd ..
cd bin2vtk/
ll
cd bin2
cd bin2vtk
ll
cd bin2vtk_xy-aho.c
cat bin2vtk_xy-aho.c
ll
cat visit_writer.c
cat bin2vtk_xy-aho.c
showq
module load git/2.9.5 
cd bec3d/
git pull
cd bec3d/
module load git/2.9.5 
git pull
cd bec3d/
module load git/2.9.5 
git add .
git commit -m "Glatka"
git pull
cat Par.py
cd MojOutput/
qsub /home/denis/bec3d/job.pbs 
ll
rm *.err
rm *.out
showq
qstat 637904
qstat -f  637904
qstat -f 637809
qstat -f 637851
qstat -f 638015
qstat -f 637757
qstat -f 637411
qstat -f 637805
qstat -f 637766
showq
qdel 637999
cd bec3d/MojOutput/
ll
qsub /home/denis/bec3d/job.pbs 
ll
showq
qstat -f 637974
showq
module load python/3.6.5 
pip --version
pip list
module purge
module load anaconda/3-5.1.0 
pip --version
python --version
pip list
pytorch --version
tensorflow --version
showq
cd bec3d/MojOutput/
ll
cat 638027.paradox.ipb.ac.rs.err
rm *.err
rm *.out
module avail
c++ -version
ls /usr/include
cd bec3d/MojOutput/
qsub /home/denis/bec3d/job.pbs 
ll
showq
ll
showq
ll
showq
ll
showq
ll
cat 638066.paradox.ipb.ac.rs.err
module load anaconda/3-5.1.0 
scipy --version
pip list
pip list
module load anaconda/3-5.1.0 
pip upgrade scipy
pip install scipy
module purge
module load python/3.6.5 
pip list
module purge
module avail
module load python/3.9.1 
pip list
module purge
module load anaconda/3-5.1.0 
pip list
pip install --upgrade scipy
pip list
module load anaconda/3-5.1.0 
conda list
showq
history
showq
pwd
ls
cd bec3d/
ls -la
cat Par.py
vim Par.py
ls -la
cd ..
pwd
cd ..
ls -la
showq
pwd
cd julija
ls -la
cd demo034
ls -la
clear
ls -lah
cd denis
ls -lah
cd bec3d/
ls -lah
ls -la
cd  MojOutput
ls -la
ls -lah
showq
uname -a
free -h
lscpu
free -h
sudo dmidecode --type memory
cd bec3d/MojOutput/
ll
rm *.err
rm *.out
ll
subq /home/denis/bec3d/job.pbs 
qsub /home/denis/bec3d/job.pbs 
showq
showq 
showq  | grep denis
showq | grep 'May 20'
ll
code -r 638067.paradox.ipb.ac.rs.err
cat [denis@paradox MojOutput]$ 
[denis@paradox MojOutput]$ code -r 638067.paradox.ipb.ac.rs.err
-bash: code: command not found
cat 638067.paradox.ipb.ac.rs.err
vim 638067.paradox.ipb.ac.rs.err
show scipy
which scipy
which numpy
pip list numpy
pip list | grep scipy
pip show scipy
python --version
cd ..
pwd
ls -la\
ls -la
pip show venv
pip list | grep venv
ls -la
pip show virtualenv
which python3
ls -la
module avail
module load python/3.6.5 
\\ls -la
ls -la
pip list
showq
qstat -f 637812
showq
python -m venv
module load python/3.6.5 
python -m venv
module purge
module load anaconda/3-5.1.0 
python -m venv
showq
ll
module avail
ls
showq
module purge
showq
module load anaconda/3-5.1.0 
pip list
ls
ll
showq
cd bec3d/
module load git/2.9.5 
git pull
cd bec3d/
cat Par.py 
cd MojOutput/
ll
rm *.err
rm *.out
ll
qsub /home/denis/bec3d/job.pbs 
ll
showq
ll
car real3d-rms.txt 
ll
cat 638465.paradox.ipb.ac.rs.err
rm *.err
rm *.out
git add .
git commit -m "malo sutra"
git pull
qsub /home/denis/bec3d/job.pbs 
ll
showq
cd bec3d/MojOutput/
ll
show
showq
cd bec3d/MojOutput/
cd plots/
ll
rm *.png
ll
cd ..
ll
showq
cat /home/denis/bec3d/Par.py 
car /home/denis/bec3d/input/real3d-input 
cat /home/denis/bec3d/input/real3d-input 
mkdir python3.11
ll
cd python3.11
curl https://www.python.org/ftp/python/3.11.3/Python-3.11.3.tar.xz
1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;1;112;112;1;0x1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c
ll
curl -O https://www.python.org/ftp/python/3.11.3/Python-3.11.3.tar.xz
tar -xf Python-3.11.3.tar.xz
ll
cd Python-3.11.3
./configure --prefix=$(pwd)
./configure --enable-optimizations
make
make install
python --version
cd ..
cd ,,
cd ..
module load Python-3.11.3
module load /home/denis/python3.11/Python-3.11.3
ll
cd python3.11/
ll
cd Python-3.11.3
ll
make install
gcc --version
module avail
module load c++/gnu/9.3.0  
make
cd env
source env/bit/activate
source /env/bit/activate
source /bit/activate
source env/bin/activate
source env/bin/activatce
cd ..
source env/bin/activate
cd env/
curl -O https://files.pythonhosted.org/packages/84/a9/2bf119f3f9cff1f376f924e39cfae18dec92a1514784046d185731301281/scipy-1.10.1.tar.gz
tar -xf scipy-<version>.tar.gz
tar -xf scipy-1.10.1.tar.gz
ll
cd scipy-1.10.1
ll
python3.11 setup.py install
cd ..
ll
rm -r scipy-1.10.1.tar.gz 
curl -o https://files.pythonhosted.org/packages/2c/d4/590ae7df5044465cc9fa2db152ae12468694d62d952b1528ecff328ef7fc/numpy-1.24.3.tar.gz
curl -0 https://files.pythonhosted.org/packages/2c/d4/590ae7df5044465cc9f
ll
module load git/2.9.5 
cd bec3d/MojOutput/
rm *.err
rm *.out
git pull
cd bec3d/MojOutput/
cat /home/denis/bec3d/Par.py 
qsub /home/denis/bec3d/job.pbs 
ll
cd plots/
ll
cd ..
ll
showq
ll
which openssl
rm -r env
ll
cd python3.11/
cd Python-3.11.3/
ll
make clean
./configure --enable-optimizations --with-ensurepip=install --with-openssl=/usr/bin/openssl
make -j 4
make clean
module list
module load c++/gnu/9.3.0 
./configure --enable-optimizations --with-ensurepip=install --with-openssl=/usr/bin/openssl
make -j 4
make clear
make clean
cd ..
ll
mkdir openssl
ll
cd openssl
curl -O https://www.openssl.org/source/openssl-3.1.0.tar.gz
tar -xf openssl-3.1.0.tar.xz
tar -xf openssl-3.1.0.tar.gz
ll
rm -r openssl-3.1.0.tar.gz 
cd openssl-3.1.0/
ll
cat IN
cat INSTALL.md 
./config --prefix=/path/to/installation/directory
./config --prefix=/home/denis/python3.11/openssl/
make
cat INSTALL.md 
ll
cd python3.11/
cd openssl/
ll
cd openssl-3.1.0/
ll
/Configure --prefix=/home/denis/python3.11/openssl/openssl-3.1.0/ --openssldir=/home/denis/python3.11/openssl/openssl-3.1.0/ssl
./Configure --prefix=/home/denis/python3.11/openssl/openssl-3.1.0/ --openssldir=/home/denis/python3.11/openssl/openssl-3.1.0/ssl
module list
module load c++/gnu/9.3.0 
./Configure --prefix=/home/denis/python3.11/openssl/openssl-3.1.0/ --openssldir=/home/denis/python3.11/openssl/openssl-3.1.0/ssl
make
ll
make
make clean
./Configure LIST
cd ..
ll
rm -r openssl/
showq
cd python3.11/
cd ..
python3.11 -m venv evn
cd python3.11/
cd Python-3.11.3/
ll
cd ..
ll
cd python3.11/
cd Python-3.11.3
module load c++/gnu/9.3.0  
make install
python --version
ll
cd ..
ll
python3 -m venv env
ll
source env/bin/activate
ll
cd python3.11/
ll
rm *.xz
ll
mkdir libr
ll
cd libr
pip install scipy
python --version
rm Python-3.11.3
cd ..
rm Python-3.11.3/
rm Python-3.11.3/
rm -r  Python-3.11.3/
ll
python --version
pip install --upgrade python
pip install --upgrade Python3.11.3
pip install --upgrade python
pip install --upgrade python==3.11.3
python --V
python --version
deactivate
ll
cd ..
ll
rm -r env
ll
python --version
cd python3.11/
curl -O https://www.python.org/ftp/python/3.11.3/Python-3.11.3.tar.xz
tar -xf Python-3.11.3.tar.xz
ll
cd Python-3.11.3
./configure --enable-optimizations --with-ensurepip=install
module list
make 
ll
make altinstall
make clean
./configure --enable-optimizations --with-ensurepip=install --prefix=$(pwd)
cat /proc/cpuinfo
make -j 4
make altinstall
ll
cd bin
ll
cd ..
cd
ll
python3.11 -m venv env
python --version
cd python3.11/
ll
rm Python-3.11.3.tar.xz
cd Python-3.11.3/
ls
ls -l python
cd ..
ls
cd Python-3.11.3/
python
ls
cd bin
ls
./python3.11 -m venv /home/denis/env
cd ..
ll
source ./env/bin/activate
cd py
cd python3.11/
cd Python-3.11.3/
ll
cd ..
ll
cd env
ll
pip3.11
pip --version
pip3.11 install scipy
cd /home/denis/env/lib/python3.11/site-packages/
ll
cd ..
ll
cd ..
ll
cd ..
ll
cd ..
ll
cd python3.11/
cd Python-3.11.3/
ll
cd bin
ll
cd ..
/home/denis/python3.11/Python-3.11.3/bin/pip3.11 install scipy
pip
pip --version
python --version
/home/denis/env/lib/python3.11/site-packages/pip install scipy
/home/denis/python3.11/Python-3.11.3/bin/pip3.11 download scipy
--user
deactivate
ll
python3 -m venv env1
sourde env1/bin/activate
ll
source env1/bin/activate
ll
cd env1
ll
cd bin
ll
cd ..
pip3.4 install scipy
cd ..
rm env1
rm -r  env1
rm -r env1
ll
cd /home/denis
ll
rm -r env1
ll
dpkg -s libssl-dev openssl
rpm -q libssl-devel openssl
module load python/3.6.5 
python3.6 -m venv env
ll
source ./env/bin/activate
python --version
pip
pip --version
cd env/
ll
pip install scipy
pip install --upgrade pip
pip list
cd ..
ll
rm -r python3.11
deactivate
module purge
source ./env/bin/activate
python --version
pip list
pip install matplotlib.pyplot
pip install matplotlib
pip list
deactivate
source ./env/bin/activate
cd bec3d/MojOutput/
ll
showq
qdel 638476
qsub /home/denis/bec3d/job.pbs 
ll
rm *.err
rm *.out
ll
source /home/denis/env/bin/activate
ll
cd env/
curl -O https://files.pythonhosted.org/packages/2c/d4/590ae7df5044465cc9fa2db152ae12468694d62d952b1528ecff328ef7fc/numpy-1.24.3.tar.gz
tar -xf numpy-1.24.3.tar.gz
ll
rm -r numpy-1.24.3.tar.gz 
pip list
cd numpy-1.24.3/
ll
cat README.md 
ll
cat INSTALL.rst 
module load c++/gnu/9.3.0 
python setup.py build
cd ..
ll
pip --version
pip install numpy
cd ..
deactivate
python3.6 -m venv env1
python3 -m venv env1
ll
source env1/bin/activate
cd env1
pip --version
pip install numpy
pip intall --upgrade pip
pip install --upgrade pip
pip --version
pip
pip install numpy
deactivate
cd ..
rm -r env1
ll
showq
cd bec3d/MojOutput/
ll
cat 638466.paradox.ipb.ac.rs.out
cd bec3d/MojOutput/
git add .
git commit -m "DRUGODRUGO"
git pull
qsub /home/denis/bec3d/Par.py 
ll
rm Par.py.e638478
rm Par.py.o638478 
ll
module avail
git pull
qsub /home/denis/bec3d/job2.py 
ll
cat 638479.paradox.ipb.ac.rs.err
rm *.err
rm *.out
mkdir drugo
ll
source /home/denis/env/bin/activate
cd drugo/
qsub /home/denis/bec3d/job2.py 
ll
cat 638480.paradox.ipb.ac.rs.err
git pull
rm *.err
rm *.out
qsub /home/denis/bec3d/job2.py 
ll
cat 638481.paradox.ipb.ac.rs.err
pip list
git pull
ll
rm *.err
rm *.out
qsub /home/denis/bec3d/job2.py 
ll
cat 638482.paradox.ipb.ac.rs.err
git pull
rm *.err
rm *.out
qsub /home/denis/bec3d/job2.py 
ll
cat 638483.paradox.ipb.ac.rs.err
cat 638483.paradox.ipb.ac.rs.out
rm *.err
rm *.out
cd ..
ll
cd drugo/
git pull
qsub /home/denis/bec3d/job2.py 
ll
cat 638484.paradox.ipb.ac.rs.err
ll
rm *.err
rm *.out
rm *.txt
cd ..
rm drugo
rm -r drugo
ll
which openssl
which openssl-devel
yum list installed openssl-devel
rpm -q --queryformat '%{INSTALLPREFIX}\n' openssl-devel
rpm -qi openssl-devel
ls /usr
ls /usr/include/openssl
ls /usr/lib/
which openssl-devel
showq
cd bec3d/MojOutput/
ll
cat 638477.paradox.ipb.ac.rs.err
rm *.err
rm *.out
cd bec3d/MojOutput/
qsub /home/denis/bec3d/job2.py 
ll
cat reald3d-rms.txt
cat real3d-rms.txt 
ll
